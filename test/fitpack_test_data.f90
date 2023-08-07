! **************************************************************************************************
!                                ____________________  ___   ________ __
!                               / ____/  _/_  __/ __ \/   | / ____/ //_/
!                              / /_   / /  / / / /_/ / /| |/ /   / ,<
!                             / __/ _/ /  / / / ____/ ___ / /___/ /| |
!                            /_/   /___/ /_/ /_/   /_/  |_\____/_/ |_|
!
!                                     A Curve Fitting Package
!
!   Refactored by Federico Perini, 10/6/2022
!   Based on the netlib library by Paul Dierckx
!
!   References :
!     - C. De Boor, "On calculating with b-splines", J Approx Theory 6 (1972) 50-62
!     - M. G. Cox, "The numerical evaluation of b-splines", J Inst Maths Applics 10 (1972) 134-149
!     - P. Dierckx, "Curve and surface fitting with splines", Monographs on numerical analysis,
!                    Oxford university press, 1993.
!
! **************************************************************************************************
module fitpack_test_data
    use fitpack_core, only: RKIND,pi
    implicit none
    private

    public :: dapasu
    public :: dapogr
    public :: dapola
    public :: daregr_x,daregr_y,daregr_z
    public :: daspgr_u,daspgr_v,daspgr_r
    public :: dasphe
    public :: dasurf_delta,dasurf_xyz

    real(RKIND), parameter :: dapola(*) = [ real(RKIND) :: &
         -0.010, -0.028, -0.002,  0.009,  0.018, -0.008, -0.180, -0.044,  0.067, -0.178, -0.035,  0.062, &
          0.010, -0.110,  0.013,  0.148,  0.087,  0.045,  0.007,  0.152,  0.030, -0.123,  0.147,  0.061, &
         -0.012, -0.249,  0.119, -0.224, -0.068,  0.111, -0.156, -0.145,  0.071,  0.246, -0.004,  0.096, &
          0.210, -0.187,  0.163,  0.187,  0.164,  0.103,  0.183,  0.236,  0.136,  0.161,  0.121,  0.061, &
         -0.237,  0.128,  0.145, -0.126,  0.183,  0.114, -0.292, -0.181,  0.148, -0.184, -0.303,  0.170, &
         -0.132, -0.335,  0.173, -0.150, -0.351,  0.173,  0.253, -0.227,  0.223,  0.256, -0.303,  0.306, &
          0.009, -0.327,  0.176,  0.274,  0.165,  0.154,  0.383,  0.016,  0.230,  0.105,  0.347,  0.187, &
          0.320,  0.064,  0.154, -0.263,  0.229,  0.228, -0.189,  0.267,  0.204, -0.292,  0.241,  0.267, &
         -0.218, -0.388,  0.233, -0.183, -0.409,  0.236, -0.361, -0.231,  0.216, -0.135, -0.392,  0.219, &
         -0.230, -0.386,  0.239,  0.403, -0.103,  0.299,  0.200, -0.347,  0.322,  0.086, -0.418,  0.314, &
          0.105,  0.403,  0.216,  0.103,  0.451,  0.271,  0.298,  0.327,  0.224,  0.482,  0.005,  0.307, &
          0.406,  0.157,  0.222, -0.122,  0.472,  0.381, -0.124,  0.420,  0.322, -0.478,  0.100,  0.361, &
         -0.435,  0.024,  0.278, -0.566, -0.124,  0.342, -0.232, -0.538,  0.334, -0.550, -0.003,  0.380, &
         -0.112, -0.512,  0.315, -0.250, -0.518,  0.308,  0.078, -0.510,  0.386,  0.137,  0.570,  0.345, &
          0.142,  0.525,  0.311,  0.125,  0.518,  0.308,  0.274,  0.494,  0.299,  0.229,  0.528,  0.284, &
          0.321,  0.433,  0.255, -0.520,  0.014,  0.360, -0.213, -0.622,  0.358, -0.668, -0.157,  0.408, &
         -0.207, -0.583,  0.329, -0.536, -0.422,  0.355, -0.441, -0.428,  0.311,  0.442,  0.447,  0.324, &
          0.239,  0.586,  0.339,  0.370,  0.502,  0.288,  0.426,  0.489,  0.311,  0.426,  0.428,  0.290, &
          0.274,  0.625,  0.347, -0.251, -0.709,  0.413, -0.445, -0.577,  0.339, -0.353, -0.674,  0.386, &
          0.668,  0.384,  0.378,  0.351,  0.704,  0.383,  0.370,  0.641,  0.359, -0.741, -0.377,  0.380, &
         -0.569, -0.633,  0.364, -0.760, -0.427,  0.383,  0.445,  0.742,  0.392,  0.675,  0.525,  0.376, &
          0.734,  0.349,  0.396,  0.442,  0.675,  0.370, -0.600, -0.740,  0.400, -0.614, -0.674,  0.384, &
          0.628,  0.708,  0.374,  0.699,  0.694,  0.395,  0.412, -0.273,  0.449, -0.226,  0.408,  0.425, &
          0.147, -0.496,  0.413,  0.345, -0.448,  0.643,  0.189, -0.556,  0.540,  0.527, -0.105,  0.411, &
         -0.130,  0.488,  0.419, -0.487,  0.143,  0.415, -0.331,  0.454,  0.623, -0.304,  0.405,  0.497, &
         -0.274,  0.431,  0.491, -0.614, -0.003,  0.410,  0.122, -0.602,  0.508,  0.091, -0.682,  0.565, &
          0.038, -0.632,  0.466,  0.621, -0.310,  0.798,  0.484, -0.418,  0.810,  0.440, -0.485,  0.844, &
          0.090,  0.634,  0.385, -0.463,  0.386,  0.709, -0.146,  0.639,  0.576, -0.654,  0.138,  0.583, &
         -0.513,  0.359,  0.749, -0.116,  0.630,  0.543, -0.003,  0.653,  0.459, -0.083,  0.674,  0.561, &
         -0.221, -0.738,  0.423, -0.765, -0.198,  0.452, -0.085, -0.736,  0.459, -0.767, -0.133,  0.442, &
          0.763, -0.100,  0.649,  0.670, -0.307,  0.863,  0.418, -0.681,  1.142,  0.515, -0.610,  1.264, &
          0.667, -0.330,  0.897,  0.713, -0.225,  0.748,  0.669, -0.362,  0.969,  0.064,  0.763,  0.493, &
          0.002,  0.789,  0.541,  0.134,  0.708,  0.430,  0.723,  0.086,  0.458,  0.708,  0.122,  0.440, &
         -0.529,  0.520,  1.098, -0.430,  0.675,  1.142, -0.731,  0.120,  0.631, -0.664,  0.392,  1.038, &
         -0.602,  0.380,  0.910, -0.467,  0.610,  1.123, -0.241,  0.695,  0.765, -0.015,  0.768,  0.544, &
         -0.242, -0.819,  0.444, -0.790, -0.259,  0.438, -0.835, -0.332,  0.418, -0.070, -0.832,  0.521, &
         -0.803, -0.124,  0.481,  0.307, -0.795,  0.985,  0.708, -0.552,  1.539,  0.138, -0.832,  0.723, &
          0.706, -0.449,  1.233,  0.612, -0.582,  1.416,  0.425, -0.686,  1.145,  0.764, -0.329,  0.997, &
          0.849, -0.041,  0.641,  0.869, -0.117,  0.717,  0.772,  0.287,  0.405,  0.287,  0.787,  0.428, &
          0.178,  0.879,  0.505,  0.848,  0.162,  0.503, -0.586,  0.573,  1.363, -0.857,  0.120,  0.717, &
         -0.102,  0.811,  0.663, -0.791,  0.233,  0.850, -0.891,  0.127,  0.752, -0.504,  0.622,  1.241, &
         -0.148,  0.858,  0.769, -0.857,  0.253,  0.927, -0.526,  0.714,  1.484, -0.150, -0.959,  0.551, &
         -0.873, -0.412,  0.419, -0.173, -0.920,  0.517, -0.956, -0.165,  0.534, -0.241, -0.913,  0.498, &
         -0.833, -0.512,  0.423, -0.821, -0.459,  0.425,  0.769, -0.531,  1.557,  0.916, -0.290,  1.039, &
          0.941, -0.046,  0.674,  0.877, -0.291,  1.009,  0.698, -0.570,  1.569,  0.412, -0.870,  1.299, &
          0.178, -0.952,  0.849,  0.342, -0.844,  1.108,  0.347, -0.934,  1.170,  0.089, -0.975,  0.740, &
          0.989,  0.014,  0.654,  0.826,  0.546,  0.400,  0.839,  0.345,  0.444,  0.294,  0.857,  0.436, &
          0.388,  0.902,  0.445,  0.328,  0.886,  0.459,  0.952,  0.147,  0.541, -0.218,  0.893,  0.880, &
         -0.079,  0.951,  0.716, -0.945,  0.157,  0.836, -0.928,  0.032,  0.652, -0.037,  0.910,  0.656, &
         -0.405,  0.913,  1.320, -0.639,  0.692,  1.756, -0.204,  0.952,  0.895, -0.896,  0.272,  0.996]

    real(RKIND), parameter :: dapogr(*) = [ real(RKIND) :: &
             0.805, 0.793, 0.796, 0.814, 0.822, 0.789, 0.813, 0.810, 0.789, 0.823, &
             0.796, 0.832, 0.786, 0.823, 0.788, 0.828, 0.831, 0.791, 0.846, 0.790, &
             0.783, 0.773, 0.776, 0.750, 0.810, 0.794, 0.797, 0.815, 0.793, 0.788, &
             0.778, 0.800, 0.814, 0.787, 0.835, 0.792, 0.760, 0.790, 0.786, 0.789, &
             0.746, 0.768, 0.733, 0.745, 0.756, 0.735, 0.737, 0.739, 0.742, 0.773, &
             0.782, 0.797, 0.809, 0.823, 0.784, 0.744, 0.734, 0.708, 0.717, 0.769, &
             0.699, 0.720, 0.699, 0.740, 0.737, 0.672, 0.699, 0.658, 0.690, 0.691, &
             0.725, 0.775, 0.798, 0.778, 0.771, 0.686, 0.732, 0.646, 0.667, 0.730, &
             0.634, 0.645, 0.653, 0.706, 0.680, 0.652, 0.639, 0.600, 0.599, 0.624, &
             0.656, 0.752, 0.755, 0.786, 0.714, 0.712, 0.644, 0.596, 0.632, 0.617, &
             0.595, 0.584, 0.620, 0.638, 0.580, 0.574, 0.539, 0.539, 0.560, 0.562, &
             0.569, 0.693, 0.721, 0.691, 0.663, 0.579, 0.587, 0.536, 0.519, 0.528, &
             0.432, 0.507, 0.498, 0.488, 0.479, 0.480, 0.431, 0.420, 0.433, 0.467, &
             0.466, 0.605, 0.638, 0.653, 0.628, 0.516, 0.447, 0.416, 0.440, 0.464, &
             0.323, 0.333, 0.350, 0.366, 0.310, 0.327, 0.317, 0.314, 0.307, 0.303, &
             0.394, 0.500, 0.627, 0.537, 0.524, 0.413, 0.298, 0.272, 0.280, 0.322, &
             0.195, 0.170, 0.191, 0.180, 0.174, 0.174, 0.165, 0.144, 0.137, 0.187, &
             0.265, 0.281, 0.413, 0.388, 0.296, 0.188, 0.161, 0.168, 0.150, 0.160, &
             0.805]

    real(RKIND), parameter :: dapasu(*) = [ real(RKIND) :: &
              0.718,  0.737,  0.837,  1.007,  1.206,  1.295,  1.294,  1.154,  0.980,  0.794,  0.701, &
              0.010,  0.004,  0.006, -0.015, -0.005, -0.003,  0.005,  0.013, -0.003, -0.008,  0.003, &
              0.930,  1.104,  1.247,  1.291,  1.242,  1.083,  0.892,  0.736,  0.705,  0.773,  0.914, &
              0.606,  0.615,  0.671,  0.797,  0.948,  1.038,  1.014,  0.936,  0.797,  0.655,  0.585, &
              0.430,  0.446,  0.529,  0.645,  0.718,  0.766,  0.726,  0.646,  0.545,  0.472,  0.418, &
              0.878,  1.044,  1.189,  1.216,  1.173,  1.027,  0.823,  0.711,  0.667,  0.736,  0.892, &
              0.254,  0.234,  0.232,  0.256,  0.306,  0.368,  0.404,  0.384,  0.373,  0.312,  0.260, &
              0.715,  0.753,  0.855,  0.994,  1.121,  1.167,  1.151,  1.055,  0.908,  0.781,  0.724, &
              0.759,  0.899,  1.007,  1.022,  0.957,  0.857,  0.728,  0.611,  0.573,  0.634,  0.763, &
             -0.232, -0.279, -0.345, -0.376, -0.388, -0.365, -0.335, -0.277, -0.234, -0.224, -0.229, &
              0.805,  0.806,  0.864,  0.942,  1.054,  1.098,  1.096,  1.033,  0.944,  0.856,  0.808, &
              0.571,  0.682,  0.717,  0.733,  0.687,  0.612,  0.512,  0.462,  0.446,  0.482,  0.565, &
             -0.715, -0.737, -0.813, -0.871, -0.883, -0.893, -0.890, -0.826, -0.760, -0.710, -0.710, &
              0.516,  0.505,  0.523,  0.566,  0.618,  0.656,  0.692,  0.665,  0.601,  0.555,  0.532, &
              0.282,  0.361,  0.408,  0.415,  0.376,  0.325,  0.260,  0.213,  0.228,  0.234,  0.294, &
             -0.905, -0.918, -0.979, -1.038, -1.077, -1.098, -1.086, -1.022, -0.975, -0.910, -0.899, &
              0.011, -0.030, -0.027, -0.048, -0.022, -0.003,  0.048,  0.035,  0.045,  0.034,  0.015, &
              0.011,  0.054,  0.088,  0.075,  0.055,  0.002, -0.070, -0.085, -0.096, -0.050, -0.002, &
             -0.691, -0.722, -0.773, -0.814, -0.876, -0.926, -0.904, -0.866, -0.802, -0.737, -0.732, &
             -0.535, -0.564, -0.615, -0.647, -0.669, -0.656, -0.632, -0.566, -0.533, -0.492, -0.532, &
             -0.308, -0.227, -0.227, -0.192, -0.249, -0.328, -0.380, -0.412, -0.416, -0.363, -0.300, &
             -0.243, -0.222, -0.247, -0.295, -0.324, -0.377, -0.391, -0.369, -0.335, -0.310, -0.235, &
             -0.794, -0.864, -0.930, -1.030, -1.105, -1.087, -1.047, -0.964, -0.854, -0.788, -0.783, &
             -0.563, -0.470, -0.420, -0.452, -0.508, -0.617, -0.696, -0.745, -0.733, -0.662, -0.560, &
              0.247,  0.294,  0.371,  0.400,  0.400,  0.373,  0.311,  0.249,  0.225,  0.211,  0.255, &
             -0.733, -0.788, -0.894, -1.031, -1.145, -1.173, -1.115, -1.006, -0.862, -0.769, -0.750, &
             -0.767, -0.636, -0.584, -0.611, -0.710, -0.840, -0.978, -1.019, -0.994, -0.905, -0.769, &
              0.580,  0.664,  0.813,  0.973,  1.017,  1.055,  0.940,  0.791,  0.655,  0.590,  0.590, &
             -0.401, -0.461, -0.551, -0.654, -0.733, -0.738, -0.715, -0.639, -0.532, -0.444, -0.421, &
             -0.896, -0.744, -0.670, -0.706, -0.836, -1.015, -1.165, -1.225, -1.199, -1.053, -0.897, &
              0.716,  0.802,  0.968,  1.154,  1.291,  1.297,  1.190,  1.035,  0.843,  0.742,  0.712, &
             -0.003, -0.011, -0.001,  0.014, -0.008, -0.006,  0.002, -0.001,  0.019, -0.010,  0.016, &
             -0.933, -0.798, -0.712, -0.759, -0.893, -1.084, -1.228, -1.303, -1.254, -1.096, -0.925, &
              0.614,  0.690,  0.817,  0.932,  1.022,  1.011,  0.934,  0.802,  0.665,  0.581,  0.582, &
              0.420,  0.461,  0.534,  0.640,  0.723,  0.742,  0.729,  0.623,  0.529,  0.446,  0.418, &
             -0.886, -0.752, -0.679, -0.726, -0.837, -1.017, -1.147, -1.229, -1.190, -1.051, -0.910, &
              0.261,  0.295,  0.359,  0.399,  0.390,  0.356,  0.331,  0.250,  0.195,  0.210,  0.252, &
              0.719,  0.771,  0.909,  1.043,  1.148,  1.173,  1.109,  1.002,  0.848,  0.758,  0.724, &
             -0.771, -0.660, -0.586, -0.627, -0.712, -0.871, -0.984, -1.035, -0.996, -0.900, -0.770, &
             -0.241, -0.238, -0.243, -0.283, -0.336, -0.381, -0.382, -0.372, -0.334, -0.296, -0.248, &
              0.794,  0.838,  0.948,  1.039,  1.090,  1.101,  1.062,  0.949,  0.853,  0.793,  0.781, &
             -0.575, -0.496, -0.407, -0.445, -0.511, -0.627, -0.706, -0.753, -0.732, -0.667, -0.558, &
             -0.722, -0.710, -0.749, -0.822, -0.878, -0.899, -0.903, -0.875, -0.786, -0.728, -0.706, &
              0.522,  0.547,  0.610,  0.656,  0.661,  0.657,  0.616,  0.580,  0.531,  0.502,  0.505, &
             -0.288, -0.239, -0.189, -0.206, -0.264, -0.308, -0.371, -0.407, -0.410, -0.349, -0.315, &
             -0.899, -0.918, -0.967, -1.047, -1.079, -1.105, -1.077, -1.034, -0.971, -0.913, -0.901, &
             -0.013,  0.034,  0.045,  0.033,  0.025,  0.014, -0.024, -0.034, -0.041, -0.029, -0.011, &
              0.020,  0.076,  0.073,  0.079,  0.053,  0.002, -0.062, -0.086, -0.081, -0.048, -0.012, &
             -0.705, -0.738, -0.812, -0.854, -0.897, -0.909, -0.866, -0.797, -0.747, -0.721, -0.720, &
             -0.514, -0.518, -0.531, -0.560, -0.621, -0.641, -0.663, -0.656, -0.597, -0.558, -0.525, &
              0.296,  0.357,  0.415,  0.430,  0.372,  0.318,  0.263,  0.213,  0.208,  0.249,  0.308, &
             -0.258, -0.284, -0.336, -0.373, -0.373, -0.371, -0.338, -0.267, -0.249, -0.236, -0.249, &
             -0.775, -0.790, -0.868, -0.947, -1.060, -1.119, -1.111, -1.040, -0.940, -0.850, -0.792, &
              0.573,  0.654,  0.727,  0.753,  0.712,  0.618,  0.508,  0.443,  0.409,  0.466,  0.565, &
              0.262,  0.214,  0.217,  0.261,  0.324,  0.360,  0.402,  0.397,  0.338,  0.318,  0.249, &
             -0.725, -0.745, -0.875, -0.995, -1.103, -1.181, -1.144, -1.042, -0.900, -0.785, -0.718, &
              0.762,  0.902,  1.001,  1.028,  0.972,  0.845,  0.734,  0.615,  0.591,  0.633,  0.763, &
              0.585,  0.594,  0.660,  0.807,  0.922,  1.010,  1.033,  0.948,  0.802,  0.671,  0.595, &
             -0.423, -0.467, -0.532, -0.636, -0.724, -0.741, -0.712, -0.639, -0.537, -0.456, -0.420, &
              0.890,  1.063,  1.186,  1.198,  1.187,  1.015,  0.853,  0.735,  0.681,  0.742,  0.894, &
              0.706,  0.718,  0.846,  1.024,  1.165,  1.290,  1.262,  1.160,  0.970,  0.828,  0.702, &
             -0.014,  0.005, -0.004,  0.006, -0.005,  0.000, -0.005, -0.007, -0.005, -0.006, -0.025, &
              0.921,  1.116,  1.255,  1.292,  1.230,  1.082,  0.883,  0.733,  0.690,  0.768,  0.909]


    ! Regrid test: x,y grid coordinates
    real(RKIND), parameter :: daregr_x(*) = [real(RKIND) :: -1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.5]
    real(RKIND), parameter :: daregr_y(*) = [real(RKIND) :: -1.5,-1.2,-0.9,-0.6,-0.3,0.0,0.3,0.6,0.9,1.2,1.5]

    ! Regrid: function values at the grid points
    real(RKIND), parameter :: daregr_z(*,*) = reshape([real(RKIND) :: &
            -0.0325, 0.0784, 0.0432, 0.0092, 0.1523, 0.0802, 0.0925,-0.0098, 0.0810,-0.0146,-0.0019, &
             0.1276, 0.0223, 0.0357, 0.1858, 0.2818, 0.1675, 0.2239, 0.1671, 0.0843, 0.0151, 0.0427, &
             0.0860, 0.1267, 0.1839, 0.3010, 0.5002, 0.4683, 0.4562, 0.2688, 0.1276, 0.1244, 0.0377, &
             0.0802, 0.1803, 0.3055, 0.4403, 0.6116, 0.7178, 0.6797, 0.5218, 0.2624, 0.1341,-0.0233, &
             0.1321, 0.2023, 0.4446, 0.7123, 0.7944, 0.9871, 0.8430, 0.6440, 0.4682, 0.1319, 0.1075, &
             0.2561, 0.1900, 0.4614, 0.7322, 0.9777, 1.0463, 0.9481, 0.6649, 0.4491, 0.2442, 0.1341, &
             0.0981, 0.2009, 0.4616, 0.5514, 0.7692, 0.9831, 0.7972, 0.5937, 0.4190, 0.1436, 0.0995, &
             0.0991, 0.1545, 0.3399, 0.4940, 0.6328, 0.7168, 0.6886, 0.3925, 0.3015, 0.1758, 0.0928, &
            -0.0197, 0.1479, 0.1225, 0.3254, 0.3847, 0.4767, 0.4324, 0.2827, 0.2287, 0.0999, 0.0785, &
             0.0032, 0.0917, 0.0246, 0.1780, 0.2394, 0.1765, 0.1642, 0.2081, 0.1049, 0.0493,-0.0502, &
             0.0101, 0.0297, 0.0468, 0.0221, 0.1074, 0.0433, 0.0626, 0.1436, 0.1092,-0.0232, 0.0132], &
             [size(daregr_x),size(daregr_y)])


   ! Spherical grid: u= latitudes, v=longitudes
   real(RKIND), parameter :: daspgr_u(*) = pi*0.05_RKIND*[2,4,6,7,8,9,10,12,14,16,18]
   real(RKIND), parameter :: daspgr_v(*) = pi*0.05_RKIND*[0,4,7,9,11,14,18,20,24,27,29,31,34,38]

   ! Spherical grid: coordinates
   real(RKIND), parameter :: daspgr_r(*) = [ real(RKIND) :: &
                                             0.400, 0.375, 0.449, 0.476, 0.428, 0.410, 0.414, &
                                             0.413, 0.413, 0.459, 0.453, 0.448, 0.425, 0.423, &
                                             0.399, 0.526, 0.666, 0.753, 0.718, 0.525, 0.458, &
                                             0.417, 0.536, 0.636, 0.768, 0.695, 0.564, 0.434, &
                                             0.389, 0.585, 1.182, 1.692, 1.349, 0.663, 0.419, &
                                             0.406, 0.611, 1.187, 1.676, 1.353, 0.697, 0.398, &
                                             0.346, 0.642, 1.605, 2.765, 1.750, 0.659, 0.397, &
                                             0.366, 0.638, 1.609, 2.723, 1.816, 0.670, 0.342, &
                                             0.359, 0.614, 1.850, 3.468, 1.821, 0.645, 0.366, &
                                             0.318, 0.613, 1.821, 3.486, 1.801, 0.629, 0.360, &
                                             0.299, 0.550, 1.645, 2.747, 1.483, 0.545, 0.309, &
                                             0.280, 0.581, 1.653, 2.771, 1.490, 0.557, 0.342, &
                                             0.302, 0.521, 1.281, 1.824, 1.088, 0.451, 0.273, &
                                             0.283, 0.542, 1.324, 1.796, 1.130, 0.430, 0.297, &
                                             0.262, 0.453, 0.806, 0.883, 0.654, 0.421, 0.298, &
                                             0.292, 0.475, 0.793, 0.880, 0.691, 0.382, 0.275, &
                                             0.335, 0.455, 0.619, 0.594, 0.524, 0.369, 0.297, &
                                             0.294, 0.439, 0.615, 0.631, 0.521, 0.372, 0.286, &
                                             0.390, 0.467, 0.560, 0.595, 0.544, 0.432, 0.343, &
                                             0.363, 0.515, 0.618, 0.587, 0.494, 0.418, 0.368, &
                                             0.533, 0.569, 0.612, 0.622, 0.615, 0.565, 0.548, &
                                             0.556, 0.604, 0.604, 0.634, 0.599, 0.534, 0.512]

   real(RKIND), parameter :: dasphe(*) = [ real(RKIND) :: &
                                       23.0,  13.0,   2.0,   1.0,  21.0,   9.0,  47.0,  22.0, &
                                       35.0,  38.0,  59.0,  10.0,  34.0,  37.0,  68.0,  18.0, &
                                       62.0,   5.0,  63.0,  20.0,  70.0,  32.0,  81.0,  28.0, &
                                      119.0,  24.0, 113.0,   4.0, 103.0,  41.0, 107.0,  19.0, &
                                      119.0,  22.0, 140.0,  18.0, 126.0,  15.0, 138.0,   2.0, &
                                      134.0,   6.0, 166.0,  29.0, 180.0,  21.0, 168.0,  34.0, &
                                       16.0,  81.0,  25.0,  87.0,  14.0,  67.0,  31.0,  61.0, &
                                       47.0,  89.0,  42.0,  72.0,  42.0,  46.0,  68.0,  60.0, &
                                       81.0,  76.0,  90.0,  84.0,  72.0,  55.0,  64.0,  57.0, &
                                      101.0,  65.0,  92.0,  54.0, 110.0,  70.0,  93.0,  90.0, &
                                      113.0,  72.0, 127.0,  53.0, 142.0,  82.0, 146.0,  49.0, &
                                      136.0,  63.0, 155.0,  75.0, 163.0,  89.0, 169.0,  87.0, &
                                        1.0, 102.0,  26.0, 120.0,  10.0,  99.0,  47.0, 102.0, &
                                       31.0, 127.0,  53.0, 120.0,  54.0, 118.0,  62.0, 100.0, &
                                       64.0, 102.0,  88.0, 106.0,  76.0, 112.0,  71.0, 129.0, &
                                      102.0, 102.0,  91.0, 132.0, 107.0, 125.0, 110.0, 103.0, &
                                      110.0, 110.0, 141.0, 133.0, 142.0, 112.0, 144.0,  93.0, &
                                      131.0, 103.0, 164.0, 118.0, 174.0,  97.0, 161.0, 109.0, &
                                       20.0, 157.0,   4.0, 158.0,  25.0, 142.0,  39.0, 164.0, &
                                       40.0, 135.0,  49.0, 177.0,  42.0, 150.0,  64.0, 147.0, &
                                       90.0, 151.0,  74.0, 168.0,  85.0, 143.0,  76.0, 158.0, &
                                      116.0, 163.0, 107.0, 160.0, 106.0, 154.0, 106.0, 152.0, &
                                       90.0, 172.0, 140.0, 163.0, 127.0, 151.0, 124.0, 144.0, &
                                      125.0, 174.0, 157.0, 140.0, 155.0, 152.0, 160.0, 165.0, &
                                       28.0, 224.0,  26.0, 220.0,  17.0, 191.0,  59.0, 196.0, &
                                       49.0, 205.0,  34.0, 216.0,  48.0, 209.0,  82.0, 219.0, &
                                       73.0, 211.0,  60.0, 203.0,  89.0, 202.0,  86.0, 191.0, &
                                       99.0, 181.0, 111.0, 214.0, 115.0, 204.0,  96.0, 185.0, &
                                       93.0, 187.0, 126.0, 188.0, 135.0, 204.0, 148.0, 190.0, &
                                      133.0, 209.0, 177.0, 207.0, 156.0, 219.0, 171.0, 199.0, &
                                       26.0, 242.0,   8.0, 263.0,  25.0, 264.0,  36.0, 239.0, &
                                       55.0, 245.0,  47.0, 225.0,  55.0, 240.0,  71.0, 232.0, &
                                       80.0, 254.0,  69.0, 236.0,  82.0, 259.0,  89.0, 252.0, &
                                      108.0, 242.0,  94.0, 229.0,  95.0, 227.0, 115.0, 250.0, &
                                      115.0, 228.0, 124.0, 236.0, 132.0, 253.0, 133.0, 263.0, &
                                      132.0, 228.0, 178.0, 267.0, 161.0, 230.0, 179.0, 234.0, &
                                       16.0, 273.0,  24.0, 296.0,   4.0, 294.0,  34.0, 309.0, &
                                       35.0, 274.0,  53.0, 279.0,  58.0, 315.0,  65.0, 296.0, &
                                       73.0, 307.0,  60.0, 276.0,  81.0, 306.0,  75.0, 277.0, &
                                      100.0, 298.0, 106.0, 290.0,  94.0, 276.0, 117.0, 293.0, &
                                       93.0, 310.0, 143.0, 296.0, 121.0, 285.0, 126.0, 282.0, &
                                      141.0, 312.0, 163.0, 306.0, 166.0, 288.0, 159.0, 302.0, &
                                       14.0, 340.0,  24.0, 353.0,  15.0, 328.0,  57.0, 333.0, &
                                       47.0, 351.0,  53.0, 330.0,  46.0, 340.0,  78.0, 326.0, &
                                       74.0, 336.0,  61.0, 339.0,  76.0, 337.0,  76.0, 356.0, &
                                      105.0, 317.0,  97.0, 343.0, 112.0, 352.0, 107.0, 316.0, &
                                      101.0, 318.0, 124.0, 340.0, 124.0, 356.0, 146.0, 346.0, &
                                      140.0, 318.0, 168.0, 327.0, 168.0, 331.0, 153.0, 333.0]

   real(RKIND), parameter :: dasurf_delta = 0.102569e-01 ! Estimate of the standard deviation
   real(RKIND), parameter :: dasurf_xyz(*,*) = reshape([ real(RKIND) :: &
                                          -1.9867, -1.9470,  0.0076, &
                                          -1.5042, -1.4471, -0.0035, &
                                          -1.7742, -1.0253,  0.0144, &
                                          -1.6004, -1.7173,  0.0036, &
                                          -1.5655, -1.8963,  0.0157, &
                                          -1.3779, -0.3124,  0.1309, &
                                          -1.9129, -0.2167,  0.0402, &
                                          -1.5747, -0.4469,  0.0625, &
                                          -1.6036, -0.3483,  0.0553, &
                                          -1.9530, -0.4068,  0.0133, &
                                          -1.3627,  0.5334,  0.1119, &
                                          -1.2517,  0.5905,  0.1269, &
                                          -1.0566,  0.8297,  0.1708, &
                                          -1.6402,  0.0373,  0.0638, &
                                          -1.3958,  0.4089,  0.1077, &
                                          -1.6710,  1.2759,  0.0274, &
                                          -1.1381,  1.9973, -0.0092, &
                                          -1.7332,  1.1393,  0.0169, &
                                          -1.7920,  1.0699,  0.0178, &
                                          -1.0052,  1.0919,  0.1106, &
                                          -0.8194, -1.7594,  0.0154, &
                                          -0.9082, -1.1300,  0.1125, &
                                          -0.3390, -1.8467,  0.0388, &
                                          -0.0594, -1.3774,  0.1521, &
                                          -0.3295, -1.1278,  0.2237, &
                                          -0.7246, -0.4482,  0.4845, &
                                          -0.3325, -0.2274,  0.8418, &
                                          -0.7718, -0.9464,  0.2268, &
                                          -0.6555, -0.0687,  0.6497, &
                                          -0.1243, -0.6415,  0.6507, &
                                          -0.0757,  0.0176,  1.0069, &
                                          -0.7134,  0.6721,  0.3690, &
                                          -0.0422,  0.6835,  0.6235, &
                                          -0.9322,  0.8158,  0.2226, &
                                          -0.3798,  0.4549,  0.7087, &
                                          -0.5879,  1.3671,  0.1145, &
                                          -0.4919,  1.1213,  0.2194, &
                                          -0.3558,  1.3006,  0.1630, &
                                          -0.0638,  1.6276,  0.0790, &
                                          -0.7318,  1.1287,  0.1610, &
                                           0.7348, -1.5365,  0.0682, &
                                           0.5143, -1.4583,  0.0947, &
                                           0.5456, -1.4447,  0.0775, &
                                           0.3855, -1.4510,  0.1061, &
                                           0.6197, -1.3444,  0.1170, &
                                           0.2048, -0.8807,  0.4397, &
                                           0.1138, -0.7657,  0.5560, &
                                           0.6618, -0.6786,  0.3930, &
                                           0.8870, -0.1306,  0.4664, &
                                           0.3598, -0.0350,  0.8872, &
                                           0.7982,  0.1373,  0.5129, &
                                           0.8622,  0.7410,  0.2846, &
                                           0.5027,  0.0037,  0.7766, &
                                           0.2402,  0.8623,  0.4501, &
                                           0.7488,  0.7138,  0.3513, &
                                           0.1189,  1.2044,  0.2338, &
                                           0.1752,  1.1834,  0.2242, &
                                           0.2245,  1.9460,  0.0087, &
                                           0.6075,  1.8895,  0.0106, &
                                           0.4312,  1.7096,  0.0472, &
                                           1.8552, -1.7396, -0.0076, &
                                           1.3081, -1.7975,  0.0080, &
                                           1.0085, -1.4331,  0.0401, &
                                           1.1658, -1.6434,  0.0157, &
                                           1.2353, -1.6872,  0.0259, &
                                           1.6036, -0.0301,  0.0826, &
                                           1.0780, -0.8749,  0.1429, &
                                           1.5563, -0.1535,  0.0760, &
                                           1.2505, -0.8546,  0.0995, &
                                           1.8409, -0.5623,  0.0319, &
                                           1.9347,  0.9197,  0.0138, &
                                           1.2037,  0.9823,  0.0934, &
                                           1.5779,  0.7897,  0.0479, &
                                           1.8671,  0.1472,  0.0276, &
                                           1.1892,  0.7828,  0.1353, &
                                           1.5763,  1.1683,  0.0123, &
                                           1.1368,  1.0035,  0.1042, &
                                           1.1464,  1.4909,  0.0203, &
                                           1.4444,  1.5232,  0.0108, &
                                           1.5092,  1.9105,  0.0022], [3,80])
end module fitpack_test_data
