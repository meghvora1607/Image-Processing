# import numpy as np
from Grid import *
from interpolate import *
from scipy.signal import convolve

def create_sinogram(phantom, num_projections, detector_sizeInPixels, detector_spacing, angular_scan_range):
    global num_proj, det_size
    det_size = detector_sizeInPixels
    num_proj = num_projections

    theta_spacing = angular_scan_range/num_projections
    sinogram = Grid(num_projections, detector_sizeInPixels, (theta_spacing, detector_spacing))
    sinogram.set_origin(0, -0.5 * (detector_sizeInPixels - 1) * detector_spacing)
    for i in range(num_projections):
        for j in range(detector_sizeInPixels):
            theta, s = sinogram.index_to_physical(i, j)
            line_integral = 0
            p = (s * np.cos(np.deg2rad(theta)), s * np.sin(np.deg2rad(theta)))
            u = (-np.sin(np.deg2rad(theta)), np.cos(np.deg2rad(theta)))
            delta_t = 0.5
            for t in np.arange(-0.5 * np.sqrt((phantom.get_size()[0]**2) + (phantom.get_size()[1]**2)) * phantom.get_spacing()[0], 0.5 * np.sqrt(2*(phantom.get_size()[0]**2) + (phantom.get_size()[1]**2)) * phantom.get_spacing()[1], 0.5):
                m = (p[0] + t * u[0], p[1] + t * u[1])
                line_integral += phantom.get_at_physical(m[0], m[1]) * delta_t
            sinogram.set_at_index(i, j, line_integral)
    return sinogram

def backproject(sinogram, reco_size_x, reco_size_y, spacing):
    reconstruction = Grid(reco_size_x, reco_size_y, spacing)
    reconstruction.set_origin(0, 0)

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

def next_power_of_two(value):
    return int(2 ** (np.ceil(np.log2(value))))


def ramp_filter(sinogram, detector_spacing):

    k = next_power_of_two(sinogram.get_size()[1])  # Length of signal after zero-padding
    delta_f = 1 / (detector_spacing * k)  # Frequency spacing
    ramp = np.abs(np.fft.fftfreq(k, delta_f))
    ramp[0:k // 2 + 1] = np.linspace(0, 1, k // 2 + 1)
    ramp[k // 2 + 1:] = np.linspace(1 - 1 / k, 0, k // 2 - 1)

    sinogram_fft = np.fft.fft(sinogram.get_buffer(), axis=1) * ramp
    sinogram_filtered = np.real(np.fft.ifft(sinogram_fft, axis=1))
    #filtered_sinogram = Grid(sinogram.get_size()[0], sinogram.get_size()[1], [detector_spacing, detector_spacing])
    sinogram.set_buffer(sinogram_filtered)
    return sinogram

def ramlak_filter(sinogram, detector_spacing):
    sinogram_buffer = sinogram.get_buffer()
    num_projections, detector_size = sinogram.get_size()
    n = detector_size // 2

    # Create the Ramachandran-Lakshminarayanan filter kernel
    kernel_size = detector_size+1
    kernel = np.zeros(kernel_size)
    j_vals = np.arange(-n, n + 1)
    kernel[n] = 0.25 / (detector_spacing ** 2)
    kernel[1::2] = -1 / (np.pi ** 2 * j_vals[1::2] ** 2 * detector_spacing ** 2)

    # Apply the filter to each projection using convolution
    filtered_sinogram = np.zeros_like(sinogram_buffer)
    for i in range(sinogram.get_size()[0]):
        projection = sinogram_buffer[i, :]
        filtered_projection = convolve(projection, kernel, mode='same')
        filtered_sinogram[i, :] = filtered_projection

    sinogram.set_buffer(filtered_sinogram)

    return sinogram