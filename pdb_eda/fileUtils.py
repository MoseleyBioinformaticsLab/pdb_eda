"""
File Utilities (pdb_eda.fileUtils)
----------------------------------

Common file utility functions used across pdb_eda's CLI.
"""

import os
import tempfile
import json

def createTempJSONFile(data, filenamePrefix):
    """Creates a temporary JSON file and returns its filename.

    :param data:  data to save into the JSON file.
    :type data: :py:class:`dict` or :py:class:`list`
    :param :py:class:`str` filenamePrefix: temporary filename prefix.
    :return: filename
    :rtype: :py:class:`str`
    """
    dirname = os.getcwd()
    filename = 0
    with tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=dirname, prefix=filenamePrefix, delete=False) as tempFile:
        json.dump(data,tempFile)
        filename = tempFile.name
    return filename

