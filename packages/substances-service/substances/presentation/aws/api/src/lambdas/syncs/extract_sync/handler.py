try:
    import os
    import sys

    # EFS - Add EFS mount path which contains libs into PYTHONPATH
    sys.path.append(os.environ["EFS_MOUNT_DIR"])  # nopep8 # noqa
except ImportError:
    pass

from substances_core.application.syncs.commands.internal.extract import (
    ExtractSyncCommand,
)

from .dependencies import extract_sync_command_handler  # type: ignore


def handle(event, context):  # type: ignore
    sync_id = event["sync_id"]

    handler = extract_sync_command_handler()
    result = handler.execute(ExtractSyncCommand(sync_id=sync_id))

    return dict(result)
