try:
    import os
    import sys

    # EFS - Add EFS mount path which contains libs into PYTHONPATH
    sys.path.append(os.environ["EFS_MOUNT_DIR"])  # nopep8 # noqa
except ImportError:
    pass

from substances_core.application.syncs.commands.internal.transform import (
    TransformSyncCommand,
)

from .dependencies import transform_sync_command_handler


def handle(event, context):
    print("Transform Sync Handler")
    sync_id = event["sync_id"]

    handler = transform_sync_command_handler()
    result = handler.execute(TransformSyncCommand(sync_id=sync_id))

    return dict(result)
