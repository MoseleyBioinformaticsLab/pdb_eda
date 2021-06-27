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
    :type data: :py:class:`dict`, :py:class:`list`
    :param filenamePrefix: temporary filename prefix.
    :type filenamePrefix: :py:class:`str`

    :return: filename
    :rtype: :py:class:`str`
    """
    dirname = os.getcwd()
    filename = 0
    with tempfile.NamedTemporaryFile(mode='w', buffering=1, dir=dirname, prefix=filenamePrefix, delete=False) as tempFile:
        json.dump(data,tempFile)
        filename = tempFile.name
    return filename

