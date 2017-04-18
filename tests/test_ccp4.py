from densityCheck import ccp4
import numpy as np
from pytest import approx

pdbid = '1cbs_diff'  # Fo - Fc
densityObj = ccp4.readFromPDBID(pdbid)
posGreen = [6.1, 21.3, 20.7]  # green, big positive
posRed = [34.3, 25.3, 22.8]  # red, big negative


def test_working():
    assert 1 + 2 == 3.0
    assert [1, 2, 3] == [1, 2, 3]
    assert [1, 2, 3] != [3, 2, 1]


def test_origin_match_LiteMol():
    """Test if the origin is the same as LiteMol calculated."""
    alpha = np.pi / 180 * densityObj.header.alpha
    beta = np.pi / 180 * densityObj.header.beta
    gamma = np.pi / 180 * densityObj.header.gamma

    xscale = densityObj.header.gridLength[0]
    yscale = densityObj.header.gridLength[1]
    zscale = densityObj.header.gridLength[2]

    z1 = np.cos(beta)
    z2 = (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    z3 = np.sqrt(1.0 - z1 * z1 - z2 * z2)

    xAxis = [xscale, 0.0, 0.0]
    yAxis = [np.cos(gamma) * yscale, np.sin(gamma) * yscale, 0.0]
    zAxis = [z1 * zscale, z2 * zscale, z3 * zscale]

    if densityObj.header.futureUse[-3] == 0.0 and densityObj.header.futureUse[-2] == 0.0 and densityObj.header.futureUse[-1] == 0.0:
        origin = [
            xAxis[0] * densityObj.header.crsStart[densityObj.header.map2xyz[0]] + yAxis[0] * densityObj.header.crsStart[densityObj.header.map2xyz[1]] +
            zAxis[0] * densityObj.header.crsStart[densityObj.header.map2xyz[2]],
            yAxis[1] * densityObj.header.crsStart[densityObj.header.map2xyz[1]] + zAxis[1] * densityObj.header.crsStart[densityObj.header.map2xyz[2]],
            zAxis[2] * densityObj.header.crsStart[densityObj.header.map2xyz[2]]
        ]
    else:
        origin = [densityObj.header.originEM[0], densityObj.header.originEM[1], densityObj.header.originEM[2]]

    assert {x == approx(y) for x, y in zip(origin, densityObj.origin)}


def test_xyz_crs_conversion():
    """Test that converting from crs to xyz and back to crs, or converting from xyz to crs back to xyz will give be the same as originals"""
    assert [60, 60, 60] == approx(densityObj.header.xyz2crsCoord(densityObj.header.crs2xyzCoord([60, 60, 60])))
    assert [80, 60, 60] == approx(densityObj.header.xyz2crsCoord(densityObj.header.crs2xyzCoord([80, 60, 60])))  # column out of boundary
    assert [60, 160, 60] == approx(densityObj.header.xyz2crsCoord(densityObj.header.crs2xyzCoord([60, 160, 60])))  # Row out of boundary
    assert [60, 60, 160] == approx(densityObj.header.xyz2crsCoord(densityObj.header.crs2xyzCoord([60, 60, 160])))  # Section out of boundary

    assert posRed == approx(densityObj.header.crs2xyzCoord(densityObj.header.xyz2crsCoord(posRed)), abs=np.max(densityObj.header.gridLength)/2)
    assert posGreen == approx(densityObj.header.crs2xyzCoord(densityObj.header.xyz2crsCoord(posGreen)), abs=np.max(densityObj.header.gridLength)/2)


def test_aberrant_point():
    """Test red and green are both beyond 3 standard deviations of overall all density"""
    assert densityObj.getPointDensityFromXyz(posRed) < - 3 * densityObj.header.rmsd
    assert densityObj.getPointDensityFromXyz(posGreen) > 3 * densityObj.header.rmsd


def test_crs_edge_cases():
    """Test some edge cases of crs coordinates, and check if they were handled correctly."""
    assert densityObj.getPointDensityFromCrs([8, 1, 30]) == densityObj.getPointDensityFromCrs([8, 77, 30])  # Data repeating after given interval number
    assert densityObj.getPointDensityFromCrs([88, 77, 30]) == densityObj.getPointDensityFromCrs([8, 77, 30])  # Out of crs boundary, so bring it back n interval number
    assert densityObj.getPointDensityFromCrs([78, 77, 30]) == 0  # If a crs coordinate was not provided in the data, set it to 0


def test_aberrant_blob():
    """Test the findAberrantBlobs funtion"""
    densityObj.density[:20, :20, :20] = 0
    centerXYZ = densityObj.header.crs2xyzCoord([10, 10, 10])
    densityCutoff = 3 * densityObj.header.rmsd

    # Set up green blobs
    densityObj.density[11:13, 11:13, 11:13] = 1  # (x,y,z), Volume 2*2*2=8
    densityObj.density[7:10, 11:14, 7:10] = 1  # (-x,y,-z), Volume 3:3:3=27

    # Set up red blobs
    densityObj.density[11:13, 7:9, 7:9] = -1  # (-x,-y,z), Volume 2*2*2=8
    densityObj.density[7:10, 7:10, 11:14] = -1  # (x,-y,-z), Volume 3*3*3=27

    resultGreen = densityObj.findAberrantBlobs(centerXYZ, 5, densityCutoff)
    centroid1 = [(densityObj.header.crs2xyzCoord([11, 11, 11])[i] + densityObj.header.crs2xyzCoord([12, 12, 12])[i])/2 for i in range(3)]
    trueGreen = [ccp4.DensityBlob(centroid1, 8, densityObj.header.unitVolume * 8, densityObj.header),
                 ccp4.DensityBlob(densityObj.header.crs2xyzCoord([8, 12, 8]), 27, densityObj.header.unitVolume * 27, densityObj.header)]

    resultGreen.sort(key=lambda x: (x.volume, x.totalDensity))
    trueGreen.sort(key=lambda x: (x.volume, x.totalDensity))
    for i in range(0, len(resultGreen)):
        assert resultGreen[i] == trueGreen[i]

    resultRed = densityObj.findAberrantBlobs(centerXYZ, 5, -1 * densityCutoff)
    centroid2 = [(densityObj.header.crs2xyzCoord([7, 7, 11])[i] + densityObj.header.crs2xyzCoord([8, 8, 12])[i]) / 2 for i in range(3)]
    trueRed = [ccp4.DensityBlob(centroid2, -8, densityObj.header.unitVolume * 8, densityObj.header),
               ccp4.DensityBlob(densityObj.header.crs2xyzCoord([12, 8, 8]), -27, densityObj.header.unitVolume * 27, densityObj.header)]

    resultRed.sort(key=lambda x: (x.volume, x.totalDensity))
    trueRed.sort(key=lambda x: (x.volume, x.totalDensity))

    for i in range(0, len(resultRed)):
        assert resultRed[i] == trueRed[i]



