from Grid import *
from Grid2 import next_power_of_two


def create_fanogram(phantom, num_projections, detector_sizeInPixels, detector_spacing, angular_increment, d_si, d_sd):
    beta_spacing = angular_increment
    fanogram = Grid(num_projections, detector_sizeInPixels, (beta_spacing, detector_spacing))
    fanogram.set_origin(0, -0.5 * (detector_sizeInPixels - 1) * detector_spacing)
    delta_t = 0.5
    for i in range(num_projections):
        for j in range(detector_sizeInPixels):
            g_t_beta = 0
            beta, t = fanogram.index_to_physical(i, j)
            source_point = -d_si*np.sin(np.deg2rad(beta)), d_si*np.cos(np.deg2rad(beta))
            detector_point = (d_sd-d_si)*np.sin(np.deg2rad(beta)), -(d_sd-d_si)*np.cos(np.deg2rad(beta))
            p = detector_point[0] + t * np.cos(np.deg2rad(beta)), detector_point[1] + t * np.sin(np.deg2rad(beta))
            sp = (p[0]-source_point[0], p[1]-source_point[1])
            length = np.sqrt(sp[1]**2 + sp[0]**2)
            no_of_samples = length/delta_t
            cosine_weighting = d_sd/(np.sqrt(d_sd**2+t**2))
            for k in range(0, int(no_of_samples)):
                sample_position = k*sp[0]/length + source_point[0], k*sp[1]/length + source_point[1]
                g_t_beta += phantom.get_at_physical(sample_position[0], sample_position[1])*delta_t*cosine_weighting
            fanogram.set_at_index(i, j, g_t_beta)
    return fanogram

def rebinning(fanogram, d_si, d_sd):
    num_proj, det_size = fanogram.get_size()
    angular_spacing, detector_spacing = fanogram.get_spacing()
    length_proj = (int)(180/angular_spacing)
    rebinning_grid = Grid((int(180/angular_spacing)), det_size, (angular_spacing, detector_spacing))
    rebinning_grid.set_origin(0, -0.5*(det_size-1)*detector_spacing)
    # value = 0
    for i in range(length_proj): #for i in range(num_proj):
        for j in range(det_size):
            theta, s = rebinning_grid.index_to_physical(i, j)
            gamma = np.arcsin(s/d_si)
            gamma_deg = np.rad2deg(gamma)
            t = d_sd * np.tan(gamma)
            beta = theta - gamma_deg
            p_val = fanogram.get_at_physical(beta, t)
            if beta<0:
                t2 = -t
                beta2 = beta + 2*gamma_deg + 180
                p_val = fanogram.get_at_physical(beta2, t2)
            rebinning_grid.set_at_index(i, j, p_val)
    return rebinning_grid

def ramp_filter(sinogram, detector_spacing):

    k = next_power_of_two(sinogram.get_size()[1])  # Length of signal after zero-padding
    delta_f = 1 / (detector_spacing * k)  # Frequency spacing
    ramp = np.abs(np.fft.fftfreq(k, delta_f))
    ramp[0:k // 2 + 1] = np.linspace(0, 1, k // 2 + 1)
    ramp[k // 2 + 1:] = np.linspace(1 - 1 / k, 0, k // 2 - 1)

    sinogram_fft = np.fft.fft(sinogram.get_buffer(), axis=1) * ramp
    sinogram_filtered = np.real(np.fft.ifft(sinogram_fft, axis=1))
    sinogram.set_buffer(sinogram_filtered)
    return sinogram

def backprojection(sinogram, reco_size_x, reco_size_y, spacing, d_si, d_sd):
    num_proj, det_size = sinogram.get_size()
    reconstruction = Grid(reco_size_x, reco_size_y, spacing)
    reconstruction.set_origin(-(reco_size_x-1)*spacing[0]/2, -(reco_size_y-1)*spacing[1]/2)

    for i in range(reco_size_x):
        for j in range(reco_size_y):
            X_i, X_j = reconstruction.index_to_physical(i, j)
            fanogram_value = 0

            for betaindex in range(num_proj):
                beta, t = sinogram.index_to_physical(betaindex, 0)
                s = -d_si * np.sin(np.deg2rad(beta)), d_si * np.cos(np.deg2rad(beta))
                SX = (X_i - s[0], X_j - s[1])
                SQ = np.dot((SX[0], SX[1]), (np.sin(np.deg2rad(beta)), -np.cos(np.deg2rad(beta))))
                U = np.linalg.norm(SQ) / d_si
                distance_ratio_alpha = d_sd/np.linalg.norm(SQ)
                SP = distance_ratio_alpha * SX[0], distance_ratio_alpha*SX[1]
                t_val = np.dot((SP[0],SP[1]),(np.cos(np.deg2rad(beta)),np.sin(np.deg2rad(beta))))
                distance_weighting = 1/(U**2)
                fanogram_value += sinogram.get_at_physical(beta, t_val)
                value = fanogram_value*distance_weighting
            value /= 2
            reconstruction.set_at_index(reco_size_x-i-1, j, value)

    return reconstruction

def backprojection_for_rebinning(sinogram, reco_size_x, reco_size_y, spacing):
    num_proj, det_size = sinogram.get_size()
    reconstruction = Grid(reco_size_x, reco_size_y, spacing)
    reconstruction.set_origin(-(reco_size_x-1)*spacing[0]/2, -(reco_size_y-1)*spacing[1]/2)

    for i in range(reco_size_x):
        for j in range(reco_size_y):
            phy_x, phy_y = reconstruction.index_to_physical(i, j)
            val = 0

            for k in range(num_proj):
                theta, s = sinogram.index_to_physical(k, 0)
                reco_val = (phy_x * np.cos(np.deg2rad(theta))) + (phy_y * np.sin(np.deg2rad(theta)))-s
                val += interpolate(sinogram, k, reco_val/spacing[1])
            reconstruction.set_at_index(reco_size_x-i-1, j, val)

    return reconstruction

