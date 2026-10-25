from pathlib import Path


def test_reference_data_exports_y_velocity_components():
    source = Path(__file__).resolve().parents[1] / "drag_function" / "reference_data.py"
    text = source.read_text(encoding="utf-8")

    assert "v1_y = v1_raw[:min_len_r, 2]" in text
    assert "v2_y = v2_raw[:min_len_r, 2]" in text
    assert text.count("'V1_y': v1_y_final") == 2
    assert text.count("'V2_y': v2_y_final") == 2


def test_reference_data_exports_signed_velocity_components():
    source = Path(__file__).resolve().parents[1] / "drag_function" / "reference_data.py"
    text = source.read_text(encoding="utf-8")

    for column, variable in {
        "V1_x": "v1_x_final",
        "V2_x": "v2_x_final",
        "V1_y": "v1_y_final",
        "V2_y": "v2_y_final",
        "V1_z": "v1_z_final",
        "V2_z": "v2_z_final",
    }.items():
        assert f"'{column}': {variable}" in text
        assert f"'{column}': np.abs({variable})" not in text
