import asyncio
import logging

import grpc

import sync_pb2
import sync_pb2_grpc


class SyncGrpcService(sync_pb2_grpc.SyncsServicer):
    # pylint: disable=no-member
    def ListSyncs(self, request, context):
        return sync_pb2.ListResponse(count=0)

    def GetSync(self, request, context):
        return sync_pb2.ItemResponse(id="Test", provider="Test")


async def __serve() -> None:
    server = grpc.aio.server()
    sync_pb2_grpc.add_SyncsServicer_to_server(SyncGrpcService(), server)
    listen_addr = "[::]:50051"
    server.add_insecure_port(listen_addr)
    logging.info("Starting server on %s", listen_addr)
    await server.start()
    await server.wait_for_termination()


def start_grpc_server():
    asyncio.run(__serve())


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    start_grpc_server()
