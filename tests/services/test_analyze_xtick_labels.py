import pytest

from two_odd_annotator.services.analyze import _two_odd_xtick_labels


def test_xtick_labels_bold_consensus_function_mathtext():
    # 10/11 are F3H => consensus at threshold 0.9
    major_char_info = {
        "2ODD13": {
            "associated_functions": ["F3H", "FNSI"],
            "associated_characterized_baits": [
                *[f"X{i}__F3H__p__1" for i in range(10)],
                "Y__FNSI__p__1",
            ],
        }
    }

    label = _two_odd_xtick_labels(["2ODD13"], major_char_info)[0]
    assert label.startswith("$(")
    assert "\\mathbf{F3H}" in label
    assert "\\mathrm{FNSI}" in label
    assert label.endswith("\\mathrm{2ODD13}$")


def test_xtick_labels_no_consensus_falls_back_to_plain_text():
    # 1/2 are F3H => no consensus with min_size=2 and threshold=0.9
    major_char_info = {
        "2ODD13": {
            "associated_functions": ["F3H", "FNSI"],
            "associated_characterized_baits": [
                "A__F3H__p__1",
                "B__FNSI__p__1",
            ],
        }
    }

    label = _two_odd_xtick_labels(["2ODD13"], major_char_info)[0]
    assert label == "(F3H,FNSI) 2ODD13"


def test_xtick_labels_escape_underscores_in_mathtext():
    major_char_info = {
        "2ODD10": {
            "associated_functions": ["F3H", "M2H_weak"],
            "associated_characterized_baits": [
                "A__M2H_weak__p__1",
                "B__M2H_weak__p__1",
            ],
        }
    }

    label = _two_odd_xtick_labels(["2ODD10"], major_char_info)[0]
    assert "\\mathbf{M2H\\_weak}" in label

