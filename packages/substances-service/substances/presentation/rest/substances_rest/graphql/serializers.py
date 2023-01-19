import graphene


class SyncModel(graphene.ObjectType):
    id = graphene.String()
    provider = graphene.String()
    state = graphene.String()
