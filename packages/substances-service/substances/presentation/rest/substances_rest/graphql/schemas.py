from collections import namedtuple

import graphene
from substances_core.application.syncs.queries.get_sync_by_id import GetSyncByIdQuery
from substances_core.domain.syncs.objects import SyncId

from .dependencies import get_sync_by_id_query_handler
from .serializers import SyncModel

SyncValueObject = namedtuple("SyncModel", ["id", "provider", "state"])


class Query(graphene.ObjectType):
    """Query"""

    sync = graphene.Field(SyncModel, id=graphene.String(required=True))
    syncs = graphene.List(SyncModel)

    def resolve_sync(self, _, id):
        result = get_sync_by_id_query_handler().handle(
            GetSyncByIdQuery(sync_id=SyncId(id))
        )
        return SyncValueObject(
            id=str(result.sync_id), provider=result.provider, state=result.state
        )

    def resolve_syncs(
        self,
        _,
    ):
        result = get_sync_by_id_query_handler().handle(
            GetSyncByIdQuery(sync_id=SyncId(id))
        )
