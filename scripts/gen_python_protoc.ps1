#python -m grpc_tools.protoc -I../../../../protos --python_out=. --pyi_out=. --grpc_python_out=. sync.proto
#python -m grpc_tools.protoc -IC:\_dev\MicroChemLabs\grpc\substances --python_out=. --pyi_out=. --grpc_python_out=. C:\_dev\MicroChemLabs\grpc\substances\sync.proto

$root_dir = (Get-Item '.').FullName
$protos_folder = "$root_dir\grpc"
$substance_protos = "$protos_folder\substances"
$substance_protos_out = "$root_dir\packages\substances-service\substances\app\app\grpc"
$reaction_protos = "$protos_folder\reactions"
$raction_protos_out = ""

function Invoke-ProtoC([string] $proto_dir, [string] $output_path) {

    $protos = Get-ChildItem -File "$proto_dir/*.proto"
    $count = $protos.Count
    if ($count -gt 0) {
        Write-Host -ForegroundColor Cyan "Processing .proto from '$proto_dir' ($count files)..."
        & python -m grpc_tools.protoc --proto_path=$proto_dir --python_out=$substance_protos_out --pyi_out=$substance_protos_out --grpc_python_out=$substance_protos_out $protos
    }
    
}

Invoke-ProtoC -proto_dir $substance_protos -output_path $substance_protos_out