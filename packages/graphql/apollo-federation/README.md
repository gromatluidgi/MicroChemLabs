# Apollo Server - Federation Supergraph

## Docker

1. Build `docker build . -t apollo-federation`
2. Run `docker run --rm -d --name vlab-apollo-federation --hostname vlab-federation --network virtual-lab -p 3000:4000 apollo-federation`


## Build

1. Generate supergrah with Rover `rover supergraph compose --config supergraph-config.yaml > ./src/supergraph.graphql`
2. or `yarn supergraph`