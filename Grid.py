import numpy as np
from interpolate import interpolate


class Grid:
    def __init__(self, height, width, spacing):
        self.origin = np.array([0, 0])
        self.height = height
        self.width = width
        self.spacing = spacing
        self.buffer = np.zeros((height, width))
        self.x = 0
        self.y = 0
        self.physical_origin = (-0.5 * (self.height - 1) * self.spacing[0], -0.5 * (self.width - 1) * self.spacing[1])

    def set_buffer(self, buffer):
        self.buffer = buffer

    def get_buffer(self):
        return self.buffer

    def get_spacing(self):
        return self.spacing

    def set_origin(self, x, y):
        self.physical_origin = (x, y)

    def get_origin(self):
        return self.origin

    def get_size(self):
        return (self.height, self.width)

    def index_to_physical(self, i, j):
        x = self.physical_origin[0] + (i * self.spacing[0])
        y = self.physical_origin[1] + (j * self.spacing[1])
        return (x, y)

    def physical_to_index(self, x, y):
        i = (x-self.physical_origin[0]) / self.spacing[0]
        j = (y-self.physical_origin[1]) / self.spacing[1]
        return (i, j)

    def set_at_index(self, i, j, val):
        self.buffer[i, j] = val

    def get_at_index(self, i, j):
        return self.buffer[i, j]

    def get_at_physical(self, i, j):
        i_x, i_y = self.physical_to_index(i, j)
        # Interpolate the value at the point
        return interpolate(self, i_x, i_y)

    # def create_sinogram(self, phantom, num_projections, detector_size, detector_spacing, angular_scan_range):
    #     # Create a grid for the sinogram
    #     sinogram = Grid(num_projections, detector_size, (0, detector_spacing))
    #     sinogram.set_origin(0, -0.5 * (detector_size - 1) * detector_spacing)
    #
    #     for i in range(num_projections):
    #         theta = i * angular_scan_range / num_projections
    #
    #         for j in range(detector_size):
    #             s_x, s_y = sinogram.index_to_physical(i, j)
    #             #x = s * np.cos(np.deg2rad(theta)) + phantom.get_size()[0] / 2
    #             #y = s * np.sin(np.deg2rad(theta)) + phantom.get_size()[1] / 2
    #             #s_x += phantom.get_size()[0]/2
    #             #s_y += phantom.get_size()[1]/2
    #             line_integral = 0
    #
    #            # for k in np.arange(-0.5, np.sqrt(phantom.get_size()[0]**2 + phantom.get_size()[1]**2), 0.5):
    #             for k in np.arange(-0.5*detector_spacing, phantom.get_size()[0], 0.5):
    #                 p = s_x * np.cos(np.deg2rad(theta)) - s_y * np.sin(np.deg2rad(theta))  #phantom.get_size()[0]/2
    #                 q = s_x * np.sin(np.deg2rad(theta)) + s_y * np.cos(np.deg2rad(theta))  #phantom.get_size()[1]/2
    #                 x_ = int(p + k * np.cos(np.deg2rad(theta)))
    #                 y_ = int(q + k * np.sin(np.deg2rad(theta)))
    #                 x_, y_ = self.physical_to_index(x_, y_)
    #                 if x_ >= 0 and x_ < phantom.get_size()[0] and y_ >= 0 and y_ < phantom.get_size()[1]:
    #                     line_integral += interpolate(phantom, x_, y_)
    #                 #x -= np.sin(np.deg2rad(theta))
    #                 #y -= np.cos(np.deg2rad(theta))
    #
    #             sinogram.set_at_index(i, j, line_integral)
    #     #np.rot90(sinogram, k=-1)
    #     return sinogram.get_buffer()

    #def create_sinogram1(self, phantom, num_projections, detector_size, detector_spacing, angular_scan_range):
        # Create a grid for the sinogram
        # sinogram = Grid(num_projections, detector_size, (0, detector_spacing))
        # sinogram.set_origin(0, -0.5 * (detector_size - 1) * detector_spacing)
        #
        # for i in range(num_projections):
        #     theta = i * angular_scan_range / num_projections
        #     rot_matrix = np.array([[np.cos(np.deg2rad(theta)), -np.sin(np.deg2rad(theta))],
        #                            [np.sin(np.deg2rad(theta)), np.cos(np.deg2rad(theta))]])
        #
        #     for j in range(detector_size):
        #         s_x, s_y = sinogram.index_to_physical(i, j)
        #         #p, q = s_x, s_y  # Initialize p and q to s_x and s_y
        #
        #         # Rotate the ray using the rotation matrix
        #         ray = rot_matrix @ np.array([[s_x], [s_y]])
        #         line_integral = 0
        #         # Calculate the line integral along the rotated ray
        #         for k in np.arange(-0.5 * detector_spacing, phantom.get_size()[0], 0.5 * detector_spacing):
        #             x_ = int(ray[0] + k * rot_matrix[0, 0])
        #             y_ = int(ray[1] + k * rot_matrix[1, 0])
        #             x_, y_ = phantom.physical_to_index(x_, y_)
        #             if x_ >= 0 and x_ < phantom.get_size()[0] and y_ >= 0 and y_ < phantom.get_size()[1]:
        #                 line_integral += interpolate(phantom, x_, y_)
        #
        #         sinogram.set_at_index(i, j, line_integral)

        #return sinogram.get_buffer()

    #def backproject(self, sinogram, reco_size_x, reco_size_y, spacing):
    #     # Set up the output image
    #     output = np.zeros((reco_size_x, reco_size_y))
    #
    #     # Set up the sinogram grid
    #     num_projections, detector_size = self.get_size()
    #     detector_spacing = spacing
    #     sinogram_origin = self.origin
    #
    #     # Loop over projection angles
    #     for i in range(num_projections):
    #         # Compute the projection angle
    #         theta = i * 180.0 / num_projections
    #
    #         # Loop over detector elements
    #         for j in range(detector_size):
    #             # Compute the physical position of the detector element
    #             x_d, y_d = self.index_to_physical(i, j)
    #             x_d -= sinogram_origin[0]
    #             y_d -= sinogram_origin[1]
    #
    #             # Compute the physical position of the corresponding image pixel
    #             x_p = x_d * np.cos(np.deg2rad(theta)) + y_d * np.sin(np.deg2rad(theta))
    #             y_p = -x_d * np.sin(np.deg2rad(theta)) + y_d * np.cos(np.deg2rad(theta))
    #             x_p += reco_size_x / 2
    #             y_p += reco_size_y / 2
    #
    #             # Convert the physical position of the image pixel to indices
    #             x_p_idx, y_p_idx = self.physical_to_index(x_p/detector_spacing, y_p/detector_spacing)
    #
    #             # Check if the pixel is within the reconstructed image
    #             if x_p_idx >= 0 and x_p_idx < reco_size_x and y_p_idx >= 0 and y_p_idx < reco_size_y:
    #                 # Add the value from the sinogram to the corresponding pixel
    #                 output[x_p_idx, y_p_idx] += sinogram[i, j]
    #
    #     return output

    # def backproject(self, sinogram, reco_size_x, reco_size_y, spacing):
    #     # Create a grid for the reconstruction
    #     reconstruction = Grid(reco_size_x, reco_size_y, (spacing, spacing))
    #     reconstruction.set_origin(-(reco_size_x - 1) / 2 * spacing, -(reco_size_y - 1) / 2 * spacing)
    #     sinogram_x , sinogram_y = self.get_size()
    #
    #     # Loop over all projection angles
    #     for i in range(sinogram_x):
    #         theta = i * spacing + sinogram_x
    #
    #         # Loop over all detector pixels
    #         for j in range(sinogram_y):
    #             s_x, s_y = sinogram.index_to_physical(i, j)
    #
    #             # Calculate the position of the pixel in the reconstruction
    #             p = s_x * np.cos(theta) + s_y * np.sin(theta)
    #             q = -s_x * np.sin(theta) + s_y * np.cos(theta)
    #             x, y = reconstruction.physical_to_index(p, q)
    #
    #             # Add the value to the corresponding pixel in the reconstruction
    #             if x >= 0 and x < reconstruction.get_size()[0] and y >= 0 and y < reconstruction.get_size()[1]:
    #                 reconstruction.set_at_index(x, y, reconstruction.get_at_index(x, y) + sinogram.get_at_index(i, j))
    #
    #     return reconstruction.get_buffer()
