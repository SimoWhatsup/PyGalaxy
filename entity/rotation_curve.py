# coding=utf-8
__author__ = 'S. Federici (DESY)'
__version__ = '0.1.0'

import numpy as np
from scipy import ndimage


class RotationCurve:
    def __init__(self, model='Bissantz2003', peculiar_velocity_of_sun='SBD2010'):
        self.model = model
        self.peculiar_velocity_of_sun = peculiar_velocity_of_sun

        if self.model == 'Bissantz2003':
            self.r_sun = 8.  # kpc
            self.v_sun = 210.  # km s-1
            self.gmax = 165.  # +165, -165(=345) deg
            self.sol = 8  # kinematically best-fitting location
        elif self.model == 'Clemens1985':
            self.r_sun = 8.5  # kpc
            self.v_sun = 220.  # km s-1
            self.gmax = 180.  # +180, -180(=360) deg
            self.sol = 8  # kinematically best-fitting location

    def get_peculiar_velocity_of_sun(self, glon, glat):
        """
        Calculate peculiar velocity of the sun based on one of the two papers
        - "DB1998":  Dehnen & Binney, 1998, MNRAS, 298, 387
        - "SBD2010": Schoenrich, Binney & Dehnen, 2010, MNRAS, 403, 1829
        :param glon:
        :param glat:
        :return:
        """
        u0 = 0.  # radial component (velocity in/out from center)
        v0 = 0.  # tangential component
        w0 = 0.  # vertical component

        if self.peculiar_velocity_of_sun == "DB1998":
            u0 = 10.0
            v0 = 5.25
            w0 = 7.17
        elif self.peculiar_velocity_of_sun == "SBD2010":
            u0 = 11.1
            v0 = 12.24
            w0 = 7.25

        vpec = u0 * np.cos(glon) * np.cos(glat) + v0 * np.sin(glon) * np.cos(glat) + w0 * np.sin(glat)
        # print glon, glat, vpec
        # print u0 * np.cos(glon) * np.cos(glat)
        # exit(0)
        return vpec * 0.005

    def compute_model(self, args):
        """
        - Bissantz2003:
        Gas-flow simulation of the inner Galaxy using smoothed particle
        hydrodynamics (SPH) and a realistic barred gravitational potential
        derived from the observed COBE DIRBE near-IR light distribution.
        (N.Bissantz et al, 2003)

        - Clemens1985:
        Massachusetts-Stony Brook Galactic plane CO survey: The Galactic disk rotation curve
        (D.P.Clemens, The Astrophysical Journal, 295:422-436, 1985)
        """
        path = args[0]
        glon = args[1]
        # glon = radians(145) for Bi
        # glon = radians(0) for Cl
        glat = args[2]
        # glat = radians(-65) for Bi
        # glat = radians(0) for Cl
        r_proj = args[3]
        dv = args[4]
        dbin = args[5]
        N = args[6]

        if self.model == 'Bissantz2003':
            # true_dis = dbin*(0.5+arange(N))
            # r_proj = true_dis*cos(glat)

            # Read in SPH model results
            # Get Bissantz's data
            file1 = path + 'testvr1.dat'
            input = open(file1, 'r')
            bissantz = input.read().split()
            n1 = 200
            vrs = np.array(bissantz).reshape(n1, n1)
            vrs = vrs.astype(np.float32)
            # free memory
            del bissantz
            # print vrs[100,98]  #-79.1474
            # print vrs[100,99]  #-56.3561
            # print vrs[100,100] #-25.6225

            R = 0.5
            xa = -10. + 0.1 * (R + np.arange(n1))
            ya = -10. + 0.1 * (R + np.arange(n1))
            rb = np.zeros((n1, n1), dtype=np.float32)
            rb = [np.sqrt(xa[i] ** 2 + ya ** 2) for i in range(n1)]
            rb = np.array(rb)
            ia = np.where(rb > 8.0)
            vrs[ia] = 0.

            # Position of sun in SPH model and
            # unit vectors (e) to GC (tangential)
            # and l = 90 (normal) direction
            x_sun_sph = 7.518  # kpc
            y_sun_sph = -2.735  # kpc
            r_sun_sph = np.sqrt(x_sun_sph ** 2 + y_sun_sph ** 2)
            ex_tan = -x_sun_sph / r_sun_sph
            ey_tan = -y_sun_sph / r_sun_sph
            ex_norm = -ey_tan
            ey_norm = ex_tan
            xha = np.zeros(3800, dtype=np.int)
            yha = np.zeros(3800, dtype=np.int)

            # Bissantz & Gerhard velocities
            xdir = ex_tan * np.cos(glon) + ex_norm * np.sin(glon)
            ydir = ey_tan * np.cos(glon) + ey_norm * np.sin(glon)
            xha = 100 + np.floor(10 * (x_sun_sph + r_proj * xdir))
            yha = 99 - np.floor(10 * (y_sun_sph + r_proj * ydir))
            ix = np.where((xha >= 0.) & (xha <= 199.) & (yha >= 0.) & (yha <= 199.))
            cx = np.size(ix[0])
            xha = xha.astype(int)
            yha = yha.astype(int)

            # Read in rotation curve
            file2 = path + 'rotcurv4.dat'
            input = open(file2, 'r')
            rotation = input.readlines()  # 4700 lines
            rotc = np.zeros((len(rotation)), dtype=np.float32)
            drot = np.zeros((len(rotation)), dtype=np.float32)
            i = 0
            for value in rotation:
                ha = float(value.split('\n')[0].split()[0])
                rotc[i] = float(value.split('\n')[0].split()[1])
                drot[i] = float(value.split('\n')[0].split()[2])
                i = i + 1

            # Calculate peculiar velocity of the sun: Equation (6)
            vpec = self.get_peculiar_velocity_of_sun(glon, glat)

            vbgr = np.zeros(N, dtype=np.float32)
            if cx > 0:
                vbgr[ix[0]] = vrs[yha[ix[0]], xha[ix[0]]]

            # Remove data holes by linear interpolation
            dmax = np.floor(2. * self.r_sun * np.cos(glon) / dbin)
            if dmax > 6:
                vba = np.zeros(vbgr.shape)
                if vbgr[0] == 0.:
                    vbgr[0] = vpec
                idx = np.where(vbgr[1:dmax] == 0.)
                cnt = np.size(idx[0])
                while cnt > 0:
                    vba[:] = 0.
                    for k in range(cnt):
                        ia = idx[0][k] + 1
                        if vbgr[ia - 1] != 0.:
                            if vbgr[ia + 1] != 0.:
                                vba[ia] = 0.5 * (vbgr[ia - 1] + vbgr[ia + 1])
                            else:
                                vba[ia] = vbgr[ia - 1]
                        else:
                            if vbgr[ia + 1] != 0.:
                                vba[ia] = vbgr[ia + 1]
                        vbgr[ia] = vba[ia]
                    idx = np.where(vbgr[1:dmax] == 0.)
                    cnt = np.size(idx[0])
            # rad = 0.01*(0.5+arange(4700))
            # plotFunc(rad,[rotc])

            # Define distance from the GC (kpc)
            radi = np.sqrt(self.r_sun ** 2 + r_proj ** 2 - 2 * self.r_sun * r_proj * np.cos(glon))  # radi.shape = (760,)

            # Radial velocity (express rotc as a function of radi)
            x = np.floor(100 * radi).astype(int)
            c = 100. * radi - x
            rot_curve = rotc[x] + c * (rotc[x + 1] - rotc[x])
            # plotFunc(radi,[rot_curve])

            # Uncorrected radial velocity: Equation (5)
            v_lsr = (self.r_sun * rot_curve / radi - self.v_sun) * np.sin(glon) * np.cos(glat)
            # plotFunc(radi,[v_lsr])

            # Extrapolate vbgr
            pmean = self.r_sun * np.cos(glon)
            if np.cos(glon) > 0.1:
                i = int(round(40. * pmean - 5.))
                vmean = np.mean(vbgr[i - 5:i + 1])
                # vmean = vbgr[i]
                vrmean = np.mean(v_lsr[i - 5:i + 1])
                # vrmean = v_lsr[i]
                # vbgr[i:] = v_lsr[i:]+(vmean-vrmean)
                vbgr[i + 1:559] = v_lsr[i + 1:559] + (vmean - vrmean)

            # Merging - transition from inner and outer Galaxy
            Rt = 9.  # kpc
            wradi = np.where(radi > Rt, 0., (Rt - radi) / 2.)
            veff = v_lsr + (vbgr - v_lsr) * np.where(wradi > 1., 1., wradi)
            # Corrected, effective velocity: Equation (7)
            fwhm = 3
            sigma = fwhm / np.sqrt(8 * np.log(2))
            veff = ndimage.gaussian_filter(veff, sigma=sigma, order=0) - vpec
            veff[-3:] = veff[-4]

            # plotFunc(radi,[veff,v_lsr,vbgr],['Veff','Vlsr','Vbgr'],position='upper right')
            # exit(0)

            # Weights from delta veff
            dveff = np.array(veff)
            dveff[-1] = np.fabs(veff[-2] - veff[-1])
            dveff[0:N - 1] = [np.fabs(veff[i + 1] - veff[i]) for i in range(N - 1)]

            # Equation (14)
            weight_veff = np.zeros(veff.shape)
            weight_veff = np.where((dveff + 1.e-8) > dv, dv, dveff + 1.e-8)

            return veff, dveff, weight_veff

        elif self.model == 'Clemens1985':
            # true_dis = dbin*(0.5+arange(N))
            # r_proj = true_dis*cos(glat)

            # Coefficients for Rsun = 8.5 kpc and Vsun = 220 kms-1
            A = [0., 3069.81, -15809.8, 43980.1, -68287.3, 54904., -17731.]
            B = [325.0912, -248.1467, 231.87099, -110.73531, 25.073006, -2.110625]
            C = [-2342.6564, 2507.60391, -1024.068760, 224.562732, -28.4080026, 2.0697271, -0.08050808, 0.00129348]
            D = 234.88

            # Conditions
            r1 = 0.09 * self.r_sun
            r2 = 0.45 * self.r_sun
            r3 = 1.6 * self.r_sun

            # Fill rotation curve
            rad = 0.01 * (0.5 + np.arange(4700))
            rotc = np.zeros(rad.shape, dtype=np.float32)
            for k in range(np.size(rad)):
                if rad[k] < r1:
                    for i in range(len(A)):
                        rotc[k] += A[i] * pow(rad[k], i)
                elif rad[k] < r2:
                    for i in range(len(B)):
                        rotc[k] += B[i] * pow(rad[k], i)
                elif rad[k] < r3:
                    for i in range(len(C)):
                        rotc[k] += C[i] * pow(rad[k], i)
                else:
                    rotc[k] += D
            # Plot the Clemens curve
            # plotFunc(rad,[rotc])

            # Define distance from GC (kpc)
            radi = np.sqrt(self.r_sun ** 2 + r_proj ** 2 - 2 * self.r_sun * r_proj * np.cos(glon))
            # print amin(radi),amax(radi),amin(rad),amax(rad)

            # Radial velocity (express rotc as a function of radi)
            x = np.floor(100. * radi).astype(int)
            c = 100. * radi - x
            rot_curve = rotc[x] + c * (rotc[x + 1] - rotc[x])
            # plotFunc(radi,[rot_curve])

            # Uncorrected radial velocity: Equation (5)
            v_lsr = (self.r_sun * rot_curve / radi - self.v_sun) * np.sin(glon) * np.cos(glat)

            # Calculate peculiar velocity of the sun: Equation (6)
            vpec = self.get_peculiar_velocity_of_sun(glon, glat)

            # Corrected, effective velocity: Equation (7)
            fwhm = 3
            sigma = fwhm / np.sqrt(8 * np.log(2))
            veff = ndimage.gaussian_filter(v_lsr, sigma=sigma, order=0) - vpec
            veff[-3:] = veff[-4]

            # plotFunc(radi,[veff,v_lsr],['Veff','Vlsr'],position='upper right')
            # exit(0)

            # Weights from delta veff
            dveff = np.array(veff)
            dveff[-1] = np.fabs(veff[-2] - veff[-1])
            dveff[0:N - 1] = [np.fabs(veff[i + 1] - veff[i]) for i in range(N - 1)]

            # Equation (14)
            weight_veff = np.zeros(veff.shape)
            weight_veff = np.where((dveff + 1.e-8) > dv, dv, dveff + 1.e-8)

            return veff, dveff, weight_veff
