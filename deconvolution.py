# coding=utf-8
__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

import numpy as np
from scipy import ndimage
from entity.rotation_curve import RotationCurve


class Deconvolution:
    def __init__(self, Tb, Tcont, Tunab, coord, vec):
        """
        Deconvolution technique - M.Pohl, P.Englmaier, and N.Bissantz
        (The Astrophysical Journal, 677:283-291, 2008)
        Limits for Bissantz2003 rotation curve: galactic longitude < |165 deg|
        """
        self.rotation_curve = RotationCurve(model=vec[11])

        # [path_curve,self.survey,self.mosaic,self.species,lat,vel,mosaic.dy,dv,utilsConf,rmin,rmax,rotcurve,maxis]
        path = vec[0]
        survey = vec[1]
        mosaic = vec[2]
        species = vec[3]
        vel = vec[5]
        dy = vec[6]
        dv = vec[7]
        C = float(vec[8]['c'])
        Ts = float(vec[8]['tspin'])
        Tcmb = float(vec[8]['tcmb'])
        rmin = vec[9]
        rmax = vec[10]
        # rotcurve = vec[11]
        maxis = vec[12]

        annuli = len(rmin)

        if maxis == 1:
            lon = vec[4]
            lat = coord
        elif maxis == 2:
            lon = coord
            lat = vec[4]

        nlon = Tb.shape[2]
        nlat = Tb.shape[1]
        nvel = Tb.shape[0]

        # free memory
        # del vec

        # Array to store results
        cubemap = np.zeros((annuli, nlat, nlon), dtype=np.float32)
        cosb = np.cos(np.radians(lat))

        # Line properties
        sigma = 0.  # velocoty dispersion (standard deviation of the distribution)
        if species == 'HI' or species == 'HI_unabsorbed':
            sigma = 4.  # hi velocity dispersion (from LAB) [km s-1]
        elif species == 'CO':
            sigma = 3.  # co velocity dispersion [km s-1]
        elif species == 'HISA':
            sigma = 2.  # hisa velocity dispersion [km s-1]

        sigma_gas = sigma
        sigma_gas_inner = 5.  # velocity dispersion inner galaxy [km s-1]

        # Define dummy line profile and its weight
        L = 20  # half-lenght covered by the line (num of v-channels)
        ivzero = int(np.floor(L / dv))  # center of the line profile (index)
        vzero = ivzero * dv  # center of the line profile (km s-1)
        iv_vec = np.arange(2 * ivzero + 1)  # vector of velocity channels (indexes) of the line profile
        v_vec = iv_vec * dv  # vector of velocity channels (km s-1) of the line profile

        # gaussian = p[0]*exp(-0.5*power( (x-p[1])/p[2],2))

        # Define dummy line profile and its weight
        line = gaussian(v_vec, [1, vzero, sigma_gas], normalized=False)
        sigma_line = sigma_gas * np.sqrt(2. * np.pi)
        lim = line / sigma_line

        # Line profile for GC region
        line_inner = gaussian(v_vec, [1, vzero, sigma_gas_inner], normalized=False)
        sigma_line_inner = sigma_gas_inner * np.sqrt(2. * np.pi)

        # Warp parameters
        # **** rx=r-10
        rx = 0.1 * np.arange(400)
        phia = (-8. + 4. * rx - 0.182 * rx ** 2) / np.exp(rx ** 2 / 400.) / 57.3
        warpb = (9. + 197. * rx - 3.1 * rx ** 2) / 1.e3
        phib = (27.13 - 4.65 * rx + 0.125 * rx ** 2) / 57.3
        warpa = (-66. + 150. * (rx - 5.) - 0.47 * (rx - 5.) ** 2) / 1.e3
        warpc = (-70. + 171. * (rx - 5.) - 5.3 * (rx - 5.) ** 2) / 1.e3
        warpa[0:53] = 0.
        warpc[0:53] = 0.

        # Physical variables
        r_sun = self.rotation_curve.r_sun  # kpc
        v_sun = self.rotation_curve.v_sun  # km s-1
        gmax = self.rotation_curve.gmax  # +165, -165(=345) deg
        sol = self.rotation_curve.sol  # kinematically best-fitting location

        z_sun = 0.015  # kpc
        gal_radius = 20.  # kpc
        gal_thick = 1.  # kpc
        dbin = 1 / gal_radius  # 0.05 kpc
        r_scale = 10.  # radial scale [kpc]

        # Cuts
        v_offset = 10.  # velocity offset [10 km/s]
        lon_inner = 20.  # inner Galaxy longitude (|l|<=20) [deg]
        residual_line = 0.5  # threshold of residual line spectrum [K km/s]
        if species == 'HISA':
            residual_line = 0.01
        amp_frac = 0.2  # percentage of the peak value [x100 %]

        # Array definition
        N = 760
        true_dis = dbin * (0.5 + np.arange(N))  # kpc

        report = open("report_" + species + ".dat", "w")

        for l in range(nlon):
            glo_deg = lon[l]
            # print "[%i/%i] longitude: %.3f deg"%(l,nlon,lon[l])
            if abs(glo_deg) <= gmax:
                glon = np.radians(glo_deg)
                dismin = np.floor(r_sun * np.fabs(np.cos(glon)) / dbin)
                radmin = np.floor(r_sun * np.fabs(np.sin(glon)) / dbin)
                # r0 = 12+radmin
                # r1 = dismin-r0
                # r2 = dismin+r0
                r0 = radmin
                # r1 = r0-4
                # r2 = r0+4
                # print glon
                # print "px =",dismin,"py =",radmin
                # print r0,r1,r2
                # print r0*dbin,r1*dbin,r2*dbin
                # exit(0)
                for b in range(nlat):
                    # print "  %i) latitude: %.3f"%(b,lat[b])
                    gla_deg = lat[b]
                    glat = np.radians(gla_deg)

                    # If 0-pixel value in the continuum map then HISA cannot be extracted
                    if species == 'HISA' and Tcont[b, l] == 0.:
                        continue

                    # Line spectrum
                    spec = np.array(nvel)
                    spec = Tb[:, b, l]
                    spec[0] = 0.
                    spec[-1] = 0.

                    # If the spectrum is empty skip this latitude (mostly for HISA)
                    if not spec.any() != 0.:
                        continue

                    if species == 'HISA':
                        # no need to smooth because the spectrum is quite clean (few channels)
                        rspec = spec
                    else:
                        rspec = ndimage.gaussian_filter(spec, sigma=sigma_gas, order=0)

                    # Find the zero-point of the spectrum
                    if 1:
                        zero_avg = 0.
                        idx_zero_avg = np.where(rspec < 0)

                        if np.size(idx_zero_avg[0]) > 0:
                            zero_avg = np.mean(rspec[idx_zero_avg])

                        spec = spec - zero_avg
                        rspec = rspec - zero_avg

                    # Define intervals and heights: Equations (4)
                    z = z_sun + true_dis * np.sin(glat)
                    r_proj = true_dis * np.cos(glat)
                    # distance from the GC
                    radi = np.sqrt(r_sun ** 2 + r_proj ** 2 - 2 * r_sun * r_proj * np.cos(glon))

                    # Get the rotation curve from a model (Bissantz2003, Clemens1985)
                    parlist = [path, glon, glat, r_proj, dv, dbin, N]
                    veff, dveff, weight_veff = self.rotation_curve.compute_model(parlist)
                    # plotFunc(r_proj, [veff,veff2], ['Bissantz 2003','Clemens 1985'], position='lower right')
                    # exit(0)

                    wco = np.fabs(dv * sum(spec))
                    wcb = wco / sigma_line
                    wco_previous = 0
                    cnt1 = 0

                    # Start deconvolution
                    while wco > residual_line:
                        ivpeak = np.argmax(rspec)
                        vpeak = vel[ivpeak]

                        if species == 'HISA' and Tcont[b, l] > Tunab[ivpeak, b, l]:
                            break

                        ivlow = ivpeak - ivzero
                        ivhigh = ivpeak + ivzero + 1

                        const = 1  # sqrt(2.*pi)
                        if ivlow >= 0:
                            iv1 = 0
                            if ivhigh < nvel:
                                iv2 = np.size(iv_vec)
                                sigma_line = sigma_gas * np.sqrt(2. * np.pi)
                                sigma_line_inner = sigma_gas_inner * np.sqrt(2. * np.pi)
                            else:
                                iv2 = np.size(iv_vec) + (nvel - ivhigh)
                                ivhigh = nvel + 1
                                sigma_line = np.fabs(dv * sum(line[iv1:iv2])) * const
                                sigma_line_inner = np.fabs(dv * sum(line_inner[iv1:iv2])) * const
                        else:
                            iv1 = -ivlow
                            iv2 = np.size(iv_vec)
                            ivlow = 0
                            sigma_line = np.fabs(dv * sum(line[iv1:iv2])) * const
                            sigma_line_inner = np.fabs(dv * sum(line_inner[iv1:iv2])) * const

                        # Find a match between gas velocity and rotation curve
                        ivgood = np.where((vpeak >= veff) & (vpeak <= (veff + dveff)))
                        cnt_ivgood = np.size(ivgood[0])

                        # linevpeak = ones(veff.shape)*vpeak
                        # plotFunc(r_proj[0:30],[veff[0:30],veff[0:30]+dveff[0:30],linevpeak[0:30]])
                        # exit(0)

                        # Standard selection of locations
                        # -------------------------------
                        # Gas with forbidden velocity is placed in distance bins with the best matching velocity
                        # except toward the inner Galaxy (|l| < 20 deg) where for a velocity offset of more than
                        # 10 km/s to the nearest allowed velocity we accept only distance bins in the GC region.
                        # delta_v_list = fabs(veff-vpeak)
                        # if((amin(delta_v_list) > v_offset) and (fabs(glo_deg) < lon_inner)):
                        #	igalactic_rad = r1+argsort(delta_v_list[r1:r2+1])

                        delta_v_list = np.fabs(veff - vpeak)
                        if ((np.amin(delta_v_list) > v_offset)):  # and (fabs(glo_deg) < lon_inner)):
                            # igalactic_rad = r1+argsort(delta_v_list[r1:r2+1])
                            # igalactic_rad = r0+argsort(delta_v_list[r1:r2+1])
                            # if r1 > r2: r1 = r2-10
                            # if r1 < 0.: r1 = 0
                            if r0 < 4: r0 = 4
                            # igalactic_rad = argsort(delta_v_list[r0-4:r0+5])#r2+10])
                            igalactic_rad = np.arange(r0 - 4, r0 + 5)
                        # velo1 = delta_v_list[r0-4:r0+5]
                        # velo2 = delta_v_list[index2[0]:index2[-1]+1]

                        # print igalactic_rad
                        # print index2
                        # print velo1
                        # print velo2
                        # exit(0)
                        # if len(igalactic_rad) == 0:
                        #	print r0,r1,r2
                        else:
                            igalactic_rad = np.argsort(delta_v_list)

                        # delta_v_list = veff-vpeak
                        # cnt_delta = size(where(delta_v_list<0)[0])

                        # veff_sign = veff/fabs(veff)
                        # vpeak_sign = vpeak/fabs(vpeak)
                        # diff_sign = where(veff_sign != vpeak_sign)
                        # cnt_diff = size(diff_sign[0])

                        # if cnt_diff > 0:
                        # print 'here2'
                        # print veff_sign,vpeak_sign
                        #	ir_sun = floor(r_sun/dbin)
                        #	igalactic_rad = argsort(delta_v_list[ir_sun-4:ir_sun+4])
                        # elif cnt_delta > 0:
                        #	igalactic_rad = argsort(delta_v_list[])


                        # The line signal is distributed among n=sol solutions with weights
                        if (cnt_ivgood == 0):
                            roots = sol
                            ilocation = np.zeros(roots, dtype=np.float32)
                            ilocation[0:roots] = igalactic_rad[0:roots]
                        else:
                            roots = cnt_ivgood + sol
                            ilocation = np.zeros(roots, dtype=np.float32)
                            ilocation[0:cnt_ivgood] = ivgood[0][0:cnt_ivgood]
                            ilocation[cnt_ivgood:roots] = igalactic_rad[0:sol]

                        # Product of three weights (velocity,kinematic,height)
                        wa = np.zeros(roots, dtype=np.float32)
                        for i, j in enumerate(ilocation):
                            # Thickness of the gas layer - equation (15)
                            # sigma_z = 1.204*((0.06-0.04*radi[j]/r_sun)+0.095*(radi[j]/r_sun)**2) # pc
                            if species == 'HI' or species == 'HI_unabsorbed' or species == 'HISA':
                                sigma_z = (100. + 30. * radi[j] / r_sun + 90. * (radi[j] / r_sun) ** 2) * 1e-3  # kpc
                            elif species == 'CO':
                                sigma_z = (50 - 50 * radi[j] / r_sun + 90 * (radi[j] / r_sun) ** 2) * 1e-3  # kpc

                            zc = 0.
                            # Warp in the outer region of the Galaxy (r >= 11 kpc)
                            if (radi[j] >= 11.):
                                sphi = r_proj[j] * np.sin(glon) / radi[j]
                                cphi = (r_proj[j] * np.cos(glon) - r_sun) / radi[j]
                                nrx = np.floor(10. * (radi[j] - 10.))
                                sphia = sphi * np.cos(phia[nrx]) - cphi * np.sin(phia[nrx])
                                sphib = 2. * sphi * cphi * np.cos(phia[nrx]) + (sphi ** 2 - cphi ** 2) * np.sin(phia[nrx])
                                # equation (16)
                                zc = (warpa[nrx] + warpb[nrx] * sphia + warpc[nrx] * sphib) * 1e-3  # kpc

                            # Weights from height above plane
                            weight_z = gaussian(z[j], [1, zc, sigma_z], normalized=False)
                            weight_z = np.where(weight_z > np.exp(-20), weight_z, np.exp(-20))
                            # Minimize the kinematically allowed but physically unlikely placing of gas
                            weight_k = gaussian(radi[j], [1, 0, r_scale], normalized=False)

                            wa[i] = weight_veff[j] * weight_k * weight_z

                        wtot = sum(wa)

                        amp = amp_frac * rspec[ivpeak]
                        amp = np.where(wcb > amp, amp, wcb)

                        # Add intensity line (= sigma_line*amp) to cubemap
                        w2 = 0.
                        for i, j in enumerate(ilocation):
                            sigma = sigma_line
                            if radi[j] < 1.:
                                w2 += wa[i] / wtot
                                sigma = sigma_line_inner
                            wga = wa[i] / wtot

                            wamp = 0.
                            deltaV = sigma  # (sigma/sqrt(2*pi))*sqrt(8*log(2))# sigma*1.2
                            if species == 'HI' or species == 'HI_unabsorbed':
                                if amp > Ts - 5: amp = Ts - 5
                                NHI = np.log((Ts - Tcmb) / (Ts - Tcmb - amp)) * Ts * deltaV * C
                                wamp = wga * NHI
                                if np.fabs(z[j]) > 1.: radi[j] = r_sun
                            elif species == 'HISA':
                                Tc = Tcont[b, l]
                                Tu = Tunab[ivpeak, b, l]
                                NHISA = get_ampHISA(rspec[ivpeak], Tc, Tu, dy, dv, r_proj[j], sigma, vec[8]) * deltaV * C
                                if NHISA == None: break
                                wamp = wga * amp_frac * NHISA
                                if np.fabs(z[j]) > 1.: radi[j] = r_sun
                            elif species == 'CO':
                                wamp = wga * amp * sigma
                                if np.fabs(z[j]) > 0.2: radi[j] = r_sun

                            for a in range(annuli):
                                if (radi[j] >= rmin[a]) and (radi[j] < rmax[a]):
                                    cubemap[a, b, l] += wamp
                                if radi[j] > 50.:
                                    print('Distance > 50. kpc! ({})'.format(radi[j]))

                        # Subtract the results from the original spectrum
                        w1 = 1. - w2
                        line1 = line[iv1:iv2] / sigma_line
                        line2 = line_inner[iv1:iv2] / sigma_line_inner
                        spec[ivlow:ivhigh] -= (w1 * line1 + w2 * line2) * amp * sigma_line

                        if species == 'HISA':
                            rspec = spec
                        else:
                            rspec = ndimage.gaussian_filter(spec, sigma=sigma_gas, order=0)

                        wco_previous = wco
                        wco = np.fabs(dv * sum(spec))
                        wcb = wco / sigma_line
                        if wco > wco_previous:
                            wco = residual_line
                            wcb = wco / sigma_line

                        cnt1 += 1
                        if cnt1 > 600:
                            string = "\ncnt1 = %i\n" % (cnt1)
                            string += "[glo,glat] = [%.4f,%.4f] - [l,b] = [%i,%i]\n" % (glo_deg, gla_deg, l, b)
                            string += "1) wco_previous = %.3f\n" % (wco_previous)
                            string += "2) wco = %.3f\n" % (wco)
                            report.write(string)

                        # print "wco = %.3f, wco_previous = %.3f"%(wco,wco_previous)
                        # if fabs(vpeak) > 250.:
                        # plotFunc(vel,[spec,rspec],['observed','gaussian filter'], position='upper right')
                        # exit(0)
                        # zcor=total(densi(370:379,gl,gb))
                        # znorm=total(densi(0:378,gl,gb))
                        # densi(370:379,gl,gb)=0.
                        # densi(0:378,gl,gb)=densi(0:378,gl,gb)*(1.+zcor/znorm)

                        # if fabs(lat[b]) < 20.:
                        #	pa = sum(cubemap[-1,b,l])
                        #	pb = sum(cubemap[:-1,b,l])
                        #	cubemap[-1,b,l] = 0.
                        #	cubemap[:-1,b,l] *= (1.+pa/pb)

                        # print cnt1,wco
                        # plotFunc(vel,[spec,rspec],['observed','gaussian filter'], position='upper right')
                        # exit(0)
        report.close()

        # radius_max = r_sun*sin(radians(10.))
        # for l in xrange(nlon):
        #	if fabs(lon[l]) < 10.:
        #		glon = radians(lon[l])
        #		for b in xrange(nlat):
        #			glat = radians(lat[b])
        #			r_proj = true_dis*cos(glat)
        #			radi = sqrt(r_sun**2+r_proj**2-2*r_sun*r_proj*cos(glon))
        # cubemap[j,:,l] += cubemap[j,:,l]
        #			j = where(radi < radius_max)
        #			print j,size(j)
        # print radi[j]
        # exit(0)

        # Corrections for latitude
        for b in xrange(nlat):
            cubemap[:, b, :] *= cosb[b]

        return cubemap
