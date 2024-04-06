import tempfile
import zipfile
from pathlib import Path

from RTModel import RTModel


def test_ps_run():
    temporary_directory = Path(tempfile.gettempdir())
    event_zip_path = Path('events/event001.zip')
    with zipfile.ZipFile(event_zip_path, 'r') as zip_file_handle:
        zip_file_handle.extractall(temporary_directory)
    rtm = RTModel(str(temporary_directory.joinpath('event001')))
    rtm.Reader()
    rtm.InitCond()
    rtm.launch_fits('PS')
    rtm.ModelSelector('PS')
