try:
    import os
    import sys

    # EFS - Add EFS mount path which contains libs into PYTHONPATH
    sys.path.append(os.environ["EFS_MOUNT_DIR"])  # nopep8 # noqa
except ImportError:
    pass

from substances_core.application.syncs.commands.internal.load import LoadSyncCommand

from .dependencies import load_sync_command_handler


def handle(event, context):
    sync_id = event["sync_id"]

    handler = load_sync_command_handler()
    result = handler.execute(LoadSyncCommand(sync_id=sync_id))

    return dict(result)
