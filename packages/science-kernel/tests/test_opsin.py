# from py4j.java_gateway import JavaGateway, java_import
# from shared.utils.path import PathUtils


# def test_gateway():
#     bin_path = PathUtils.get_parent_dir(__file__, 1)
#     gateway = JavaGateway.launch_gateway(
#         classpath=PathUtils.join_path(
#             bin_path,
#             "libs",
#             "opsin-cli-2.7.0-jar-with-dependencies.jar",
#         ),
#         die_on_exit=True,
#     )
#     java_import(gateway.jvm, "uk.ac.cam.ch.wwmm.opsin.*")

#     svc = gateway.jvm.NameToInchi()

#     inchi = svc.parseToStdInchi("acetamide")

#     assert inchi is not None
