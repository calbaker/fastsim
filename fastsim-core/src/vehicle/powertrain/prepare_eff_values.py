"""
Python file demonstrating how to take altrios formatted efficiency arrays and
turn them into arrays that can be accepted into FASTSim-3 Interpolator with x
being power, y being SOC, and z being temperature.
"""

import numpy as np

# efficiency array as a function of power and SOC at constant 23 degrees C,
# corresponds to eta_interp_values[0] in Altrios
arr_2d = np.array([
    [
        0.760718703139553,
        0.859657826199026,
        0.953837950614283,
        0.976952217143662,
        0.995392686036557,
        0.996199812842203,
        0.980396055723986,
        0.959140335181896,
        0.9104576054872,
        0.850496066376365,
        0.662822531196789,
    ],
    [
        0.819665292388297,
        0.893333922110668,
        0.964734030200656,
        0.982382195786577,
        0.996477451881138,
        0.996869317429925,
        0.983940489578339,
        0.966787357566917,
        0.928569818352618,
        0.883690732054488,
        0.761323368524977,
    ],
    [
        0.861126295590167,
        0.917434472928836,
        0.972618799674411,
        0.986316630969094,
        0.997263804480139,
        0.997365818611071,
        0.986543220595331,
        0.972326035966636,
        0.941254620135046,
        0.905898343041854,
        0.816480936618867,
    ],
    [
        0.878691855344326,
        0.927738797288824,
        0.976009157588824,
        0.988009479867646,
        0.997602218913443,
        0.997647289208985,
        0.988009132699866,
        0.975417783858571,
        0.948187564941742,
        0.917712482642518,
        0.843364377303671,
    ],
    [
        0.886131975996052,
        0.932118886292756,
        0.977453428587132,
        0.988730795697105,
        0.997746427628651,
        0.997784899559042,
        0.988723312765708,
        0.976916942970075,
        0.951512714431768,
        0.923302706871554,
        0.855586635744082,
    ],
    [
        0.892085814015878,
        0.935630351894604,
        0.978612548997097,
        0.989309767361546,
        0.997862182940346,
        0.997883658588665,
        0.989234853052088,
        0.977987897893482,
        0.953873809588552,
        0.927243266086586,
        0.864025853742296,
    ],
    [
        0.896832630267769,
        0.938433873035159,
        0.979538756999652,
        0.989772443355399,
        0.997954689805819,
        0.997972725386608,
        0.989695472107132,
        0.978950231777395,
        0.95598541366672,
        0.930747528665383,
        0.871414333266251,
    ],
    [
        0.893575286894762,
        0.936509679799859,
        0.978902982272163,
        0.989454845828703,
        0.997891189482879,
        0.99800142651098,
        0.9898437584036,
        0.979259630067408,
        0.956662315949007,
        0.931866941724241,
        0.873752195303617,
    ],
    [
        0.901255647366545,
        0.941049211458201,
        0.980403395053744,
        0.990204395307691,
        0.998041055929327,
        0.998188861207792,
        0.990810422644317,
        0.981271771420741,
        0.961041093661913,
        0.93906309562653,
        0.888534008899869,
    ],
    [
        0.910205753769407,
        0.946350126361949,
        0.982157591234071,
        0.99108083995559,
        0.99821630174479,
        0.998303554874889,
        0.991400462700334,
        0.982495895783204,
        0.963685443586201,
        0.943371735300725,
        0.897188956161066,
    ],
    [
        0.919791284508856,
        0.952039657733531,
        0.98404276652578,
        0.992022851943828,
        0.99840466626583,
        0.998358081465211,
        0.991680583107941,
        0.983075975809841,
        0.964933439109322,
        0.945395684731206,
        0.901206383889721,
    ],
])

# transposing the altrios array so that the outermost layer is now power, and
# the innermost layer SOC (in altrios, the outermost layer is SOC and innermost
# is power), as required in order to use this array in a FASTSim-3 interpolator,
# with x being power. arr_2d_transposed can now be used in FASTSim-3.
arr_2d_transposed = np.transpose(arr_2d, (1,0))
# %%

# efficiency array as a function of power SOC, and temperature, corresponds to
# eta_interp_values in Altrios
arr_3d = [
    [
        [
            0.760718703139553,
            0.859657826199026,
            0.953837950614283,
            0.976952217143662,
            0.995392686036557,
            0.996199812842203,
            0.980396055723986,
            0.959140335181896,
            0.9104576054872,
            0.850496066376365,
            0.662822531196789,
        ],
        [
            0.819665292388297,
            0.893333922110668,
            0.964734030200656,
            0.982382195786577,
            0.996477451881138,
            0.996869317429925,
            0.983940489578339,
            0.966787357566917,
            0.928569818352618,
            0.883690732054488,
            0.761323368524977,
        ],
        [
            0.861126295590167,
            0.917434472928836,
            0.972618799674411,
            0.986316630969094,
            0.997263804480139,
            0.997365818611071,
            0.986543220595331,
            0.972326035966636,
            0.941254620135046,
            0.905898343041854,
            0.816480936618867,
        ],
        [
            0.878691855344326,
            0.927738797288824,
            0.976009157588824,
            0.988009479867646,
            0.997602218913443,
            0.997647289208985,
            0.988009132699866,
            0.975417783858571,
            0.948187564941742,
            0.917712482642518,
            0.843364377303671,
        ],
        [
            0.886131975996052,
            0.932118886292756,
            0.977453428587132,
            0.988730795697105,
            0.997746427628651,
            0.997784899559042,
            0.988723312765708,
            0.976916942970075,
            0.951512714431768,
            0.923302706871554,
            0.855586635744082,
        ],
        [
            0.892085814015878,
            0.935630351894604,
            0.978612548997097,
            0.989309767361546,
            0.997862182940346,
            0.997883658588665,
            0.989234853052088,
            0.977987897893482,
            0.953873809588552,
            0.927243266086586,
            0.864025853742296,
        ],
        [
            0.896832630267769,
            0.938433873035159,
            0.979538756999652,
            0.989772443355399,
            0.997954689805819,
            0.997972725386608,
            0.989695472107132,
            0.978950231777395,
            0.95598541366672,
            0.930747528665383,
            0.871414333266251,
        ],
        [
            0.893575286894762,
            0.936509679799859,
            0.978902982272163,
            0.989454845828703,
            0.997891189482879,
            0.99800142651098,
            0.9898437584036,
            0.979259630067408,
            0.956662315949007,
            0.931866941724241,
            0.873752195303617,
        ],
        [
            0.901255647366545,
            0.941049211458201,
            0.980403395053744,
            0.990204395307691,
            0.998041055929327,
            0.998188861207792,
            0.990810422644317,
            0.981271771420741,
            0.961041093661913,
            0.93906309562653,
            0.888534008899869,
        ],
        [
            0.910205753769407,
            0.946350126361949,
            0.982157591234071,
            0.99108083995559,
            0.99821630174479,
            0.998303554874889,
            0.991400462700334,
            0.982495895783204,
            0.963685443586201,
            0.943371735300725,
            0.897188956161066,
        ],
        [
            0.919791284508856,
            0.952039657733531,
            0.98404276652578,
            0.992022851943828,
            0.99840466626583,
            0.998358081465211,
            0.991680583107941,
            0.983075975809841,
            0.964933439109322,
            0.945395684731206,
            0.901206383889721,
        ],
    ],
    [
        [
            0.834097, 0.901686, 0.967459, 0.983741, 0.996749, 0.997053, 0.984908, 0.968854,
            0.933344, 0.892143, 0.783115,
        ],
        [
            0.865745, 0.920139, 0.973507, 0.98676, 0.997352, 0.997505, 0.987267, 0.973855,
            0.944696, 0.911791, 0.830076,
        ],
        [
            0.888862, 0.933728, 0.977985, 0.988996, 0.997799, 0.997848, 0.989051, 0.977603,
            0.953027, 0.925833, 0.861023,
        ],
        [
            0.901761, 0.941348, 0.980502, 0.990254, 0.998051, 0.998067, 0.990182, 0.979964,
            0.958199, 0.934401, 0.879006,
        ],
        [
            0.906979, 0.944437, 0.981524, 0.990764, 0.998153, 0.998183, 0.990781, 0.981211,
            0.96091, 0.938849, 0.888101,
        ],
        [
            0.91092, 0.946773, 0.982298, 0.991151, 0.99823, 0.998247, 0.991108, 0.98189, 0.962378,
            0.941245, 0.892934,
        ],
        [
            0.91391, 0.948548, 0.982885, 0.991445, 0.998289, 0.998316, 0.991465, 0.982629,
            0.963972, 0.943838, 0.898117,
        ],
        [
            0.912025, 0.947429, 0.982515, 0.991259, 0.998252, 0.998344, 0.991608, 0.982926,
            0.964611, 0.944874, 0.900173,
        ],
        [
            0.918509, 0.951278, 0.98379, 0.991897, 0.998379, 0.998484, 0.992325, 0.984407,
            0.967786, 0.949999, 0.910233,
        ],
        [
            0.925533, 0.955453, 0.985175, 0.992589, 0.998518, 0.998583, 0.992835, 0.985459,
            0.970027, 0.953595, 0.917182,
        ],
        [
            0.932697, 0.959719, 0.986591, 0.993296, 0.998659, 0.99863, 0.993074, 0.985952,
            0.971073, 0.955267, 0.920382,
        ],
    ],
    [
        [
            0.886107, 0.932104, 0.977449, 0.988728, 0.997746, 0.997831, 0.98896, 0.977413,
            0.952609, 0.925135, 0.859528,
        ],
        [
            0.903525, 0.942392, 0.980848, 0.990426, 0.998085, 0.998128, 0.990498, 0.980623,
            0.959634, 0.936759, 0.883848,
        ],
        [
            0.916625, 0.950159, 0.983419, 0.991711, 0.998342, 0.998354, 0.991659, 0.983031,
            0.964837, 0.945239, 0.900897,
        ],
        [
            0.924612, 0.954905, 0.984993, 0.992498, 0.9985, 0.99849, 0.992359, 0.984479, 0.967938,
            0.950244, 0.91071,
        ],
        [
            0.926955, 0.956299, 0.985456, 0.992729, 0.998546, 0.998581, 0.992823, 0.985434,
            0.969973, 0.953508, 0.917016,
        ],
        [
            0.929788, 0.957986, 0.986015, 0.993009, 0.998602, 0.998645, 0.993153, 0.986113,
            0.971415, 0.955812, 0.921423,
        ],
        [
            0.931874, 0.959228, 0.986428, 0.993215, 0.998643, 0.998659, 0.993222, 0.986257,
            0.97172, 0.956297, 0.922346,
        ],
        [
            0.930659, 0.958505, 0.986188, 0.993095, 0.998619, 0.998694, 0.993403, 0.986628,
            0.972504, 0.957547, 0.924718,
        ],
        [
            0.935654, 0.961481, 0.987176, 0.993589, 0.998718, 0.998778, 0.99383, 0.987505,
            0.974356, 0.960486, 0.930258,
        ],
        [
            0.941492, 0.964962, 0.988333, 0.994167, 0.998833, 0.998855, 0.994223, 0.98831,
            0.976049, 0.963164, 0.935258,
        ],
        [
            0.947771, 0.968711, 0.989579, 0.99479, 0.998958, 0.998932, 0.994615, 0.989112,
            0.977731, 0.965814, 0.940161,
        ],
    ],
    [
        [
            0.884085, 0.930913, 0.977056, 0.988532, 0.997707, 0.997719, 0.988383, 0.976204,
            0.949933, 0.920654, 0.849832,
        ],
        [
            0.901613, 0.94126, 0.980473, 0.990239, 0.998048, 0.998078, 0.990239, 0.980084,
            0.958461, 0.934833, 0.879896,
        ],
        [
            0.914783, 0.949066, 0.983057, 0.99153, 0.998306, 0.998344, 0.991609, 0.982928,
            0.964615, 0.94488, 0.900185,
        ],
        [
            0.922833, 0.953848, 0.984642, 0.992322, 0.998465, 0.998481, 0.992311, 0.984379,
            0.967726, 0.949903, 0.910046,
        ],
        [
            0.926256, 0.955884, 0.985318, 0.99266, 0.998532, 0.998535, 0.992586, 0.984947,
            0.968936, 0.951846, 0.913814,
        ],
        [
            0.929788, 0.957986, 0.986015, 0.993009, 0.998602, 0.998618, 0.993014, 0.985827,
            0.970809, 0.954845, 0.919576,
        ],
        [
            0.930879, 0.958635, 0.986231, 0.993116, 0.998623, 0.998633, 0.993088, 0.98598,
            0.971132, 0.95536, 0.920561,
        ],
        [
            0.930659, 0.958505, 0.986188, 0.993095, 0.998619, 0.998678, 0.993318, 0.986452,
            0.972134, 0.956958, 0.923601,
        ],
        [
            0.934754, 0.960944, 0.986998, 0.9935, 0.9987, 0.998786, 0.993871, 0.987588, 0.974531,
            0.960764, 0.930778,
        ],
        [
            0.940639, 0.964453, 0.988163, 0.994082, 0.998817, 0.998848, 0.994185, 0.988231,
            0.975884, 0.962903, 0.934772,
        ],
        [
            0.947127, 0.968326, 0.989451, 0.994726, 0.998945, 0.998911, 0.994506, 0.98889,
            0.977265, 0.96508, 0.938808,
        ],
    ],
]

# transposing the altrios array so that the outermost layer is now power, and
# the innermost layer temperature (in altrios, the outermost layer is
# temperature and innermost is power), as required in order to use this array in
# a FASTSim-3 interpolator, with x being power. arr_2d_transposed can now be
# used in FASTSim-3.
arr_3d_transposed = np.transpose(arr_3d, (2,1,0))