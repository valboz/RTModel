from pathlib import Path

from RTModel import RTModel


def test_reader_executable_exists():
    rt_model = RTModel()
    reader_executable_path = Path(rt_model.bindir + rt_model.readerexe)
    assert reader_executable_path.exists()


def test_init_cond_executable_exists():
    rt_model = RTModel()
    init_cond_executable_path = Path(rt_model.bindir + rt_model.initcondexe)
    assert init_cond_executable_path.exists()


def test_lev_mar_executable_exists():
    rt_model = RTModel()
    lev_mar_executable_path = Path(rt_model.bindir + rt_model.levmarexe)
    assert lev_mar_executable_path.exists()


def test_model_selector_executable_exists():
    rt_model = RTModel()
    model_selector_executable_path = Path(rt_model.bindir + rt_model.modelselectorexe)
    assert model_selector_executable_path.exists()


def test_finalizer_executable_exists():
    rt_model = RTModel()
    finalizer_executable_path = Path(rt_model.bindir + rt_model.finalizerexe)
    assert finalizer_executable_path.exists()
