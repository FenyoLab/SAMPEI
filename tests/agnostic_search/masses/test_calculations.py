from src.agnostic_search.masses.calculations import select_fragments_scan
import numpy as np


def test_select_fragment_scan():
    frag_mz_min = 200
    fragment_error = 20.0
    fragment_error_type = "ppm"
    num_intense = 20
    pep_len = 13
    precursor_error = 300.0
    precursor_error_type = "ppm"
    fragments = [
        (-127.01749984343795, "b1"),
        (-145.02806454343795, "b1o"),
        (-144.04404894843793, "b1*"),
        (-39.985471408437945, "b2"),
        (-57.99603610843795, "b2o"),
        (-57.012020513437946, "b2*"),
        (120.04518309656206, "b3"),
        (102.03461839656205, "b3o"),
        (103.01863399156206, "b3*"),
        (221.09286160156205, "b4"),
        (203.08229690156205, "b4o"),
        (204.06631249656203, "b4*"),
        (322.14054010656207, "b5"),
        (304.1299754065621, "b5o"),
        (305.11399100156206, "b5*"),
        (436.18346757656207, "b6"),
        (418.1729028765621, "b6o"),
        (419.15691847156205, "b6*"),
        (596.2141220815621, "b7"),
        (578.2035573815621, "b7o"),
        (579.1875729765621, "b7*"),
        (709.2981860965621, "b8"),
        (691.287621396562, "b8o"),
        (692.2716369915621, "b8*"),
        (780.335299901562, "b9"),
        (762.324735201562, "b9o"),
        (763.3087507965621, "b9*"),
        (877.388063776562, "b10"),
        (859.377499076562, "b10o"),
        (860.3615146715621, "b10*"),
        (990.472127791562, "b11"),
        (972.461563091562, "b11o"),
        (973.445578686562, "b11*"),
        (1061.509241596562, "b12"),
        (1043.498676896562, "b12o"),
        (1044.4826924915621, "b12*"),
        (147.11280374999998, "y1"),
        (129.10223904999998, "y1o"),
        (130.086254645, "y1*"),
        (218.149917555, "y2"),
        (200.139352855, "y2o"),
        (201.12336845, "y2*"),
        (331.23398156999997, "y3"),
        (313.22341687, "y3o"),
        (314.20743246499995, "y3*"),
        (428.28674544499995, "y4"),
        (410.276180745, "y4o"),
        (411.26019633999994, "y4*"),
        (499.32385924999994, "y5"),
        (481.31329454999997, "y5o"),
        (482.2973101449999, "y5*"),
        (612.4079232649999, "y6"),
        (594.3973585649999, "y6o"),
        (595.38137416, "y6*"),
        (772.4385777699999, "y7"),
        (754.4280130699999, "y7o"),
        (755.412028665, "y7*"),
        (886.4815052399999, "y8"),
        (868.4709405399999, "y8o"),
        (869.454956135, "y8*"),
        (987.529183745, "y9"),
        (969.5186190449999, "y9o"),
        (970.50263464, "y9*"),
        (1088.57686225, "y10"),
        (1070.56629755, "y10o"),
        (1071.550313145, "y10*"),
        (1248.607516755, "y11"),
        (1230.596952055, "y11o"),
        (1231.58096765, "y11*"),
        (1335.63954519, "y12"),
        (1317.62898049, "y12o"),
        (1318.612996085, "y12*"),
        (604.3110226732811, "[M+H]+2"),
        (595.305740323281, "[M+H-H2O]+2"),
        (595.7977481207811, "[M+H-NH3]+2"),
        (586.300457973281, "[M+H-2*H2O]+2"),
        (587.2844735682811, "[M+H-2*NH3]+2"),
        (1207.614769346562, "[M+H]+1"),
        (1189.604204646562, "[M+H-H2O]+1"),
        (1190.588220241562, "[M+H-NH3]+1"),
        (1171.593639946562, "[M+H-2*H2O]+1"),
        (1173.5616711365622, "[M+H-2*NH3]+1"),
    ]
    ions = [
        (173.128067, 28878.576171875),
        (606.359436, 24980.537109375),
        (199.1799927, 24198.08984375),
        (303.1766357, 23738.78515625),
        (120.0804749, 19147.056640625),
        (1007.497131, 18925.765625),
        (201.1229095, 16891.4765625),
        (326.1698914, 15669.6826171875),
        (772.4285278, 14136.2001953125),
        (147.1123657, 13019.501953125),
        (921.4407349, 12882.1640625),
        (175.1184235, 12809.970703125),
        (129.1019135, 12705.66796875),
        (227.1748047, 12238.0654296875),
        (602.3204346, 11504.740234375),
        (213.0864716, 10814.30859375),
        (143.1175079, 10454.1171875),
        (136.0753326, 9952.8671875),
        (701.3925781, 9830.0927734375),
        (1008.497803, 9630.6494140625),
    ]

    (
        matched_intensity,
        sequence_evidence_b,
        sequence_evidence_y,
        mgf_pepmatch_ions,
    ) = select_fragments_scan(
        frag_mz_min,
        fragments,
        precursor_error,
        precursor_error_type,
        fragment_error,
        fragment_error_type,
        num_intense,
        ions,
        pep_len,
    )

    matched_intensity_expected = 31027.6767578125
    sequence_evidence_b_expected = np.array([0.0 for x in range(13)])
    sequence_evidence_y_expected = np.array(
        [
            25725.16992188,
            16891.4765625,
            0.0,
            0.0,
            0.0,
            0.0,
            14136.20019531,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    mgf_pepmass_ions_expected = [
        ["y1", 147.11280374999998, 147.1123657, 13019.501953125],
        ["y1o", 129.10223904999998, 129.1019135, 12705.66796875],
        ["y2*", 201.12336845, 201.1229095, 16891.4765625],
        ["y7", 772.4385777699999, 772.4285278, 14136.2001953125],
    ]
    assert matched_intensity == matched_intensity_expected
    assert all(
        [
            abs(a - b) < 1e-4
            for a, b in zip(sequence_evidence_b, sequence_evidence_b_expected)
        ]
    )
    assert all(
        [
            abs(a - b) < 1e-4
            for a, b in zip(sequence_evidence_y, sequence_evidence_y_expected)
        ]
    )
    assert all(
        [
            all([a == b for a, b in zip(returned, expected)])
            for returned, expected in zip(mgf_pepmatch_ions, mgf_pepmass_ions_expected)
        ]
    )


# matched_intensity = 31027.6767578125
# matched_intensity_parent = 0.0
# matched_intensity_not_parent = 31027.6767578125
