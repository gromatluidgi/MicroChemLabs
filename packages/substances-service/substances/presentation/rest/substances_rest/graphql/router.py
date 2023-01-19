import graphene
from fastapi import APIRouter
from starlette_graphene3 import GraphQLApp, make_graphiql_handler

from .schemas import Query

graphql_router = APIRouter()

graphql_router.add_route(
    "/graphql",
    GraphQLApp(schema=graphene.Schema(query=Query), on_get=make_graphiql_handler()),
)
