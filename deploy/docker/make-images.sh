# Makes the base dependencies
docker build --no-cache -f 001b-mpich-3.2-yandabase -t yandabase .
docker build --no-cache -f 001b-depends -t yandabase .
# Makes yandasoft
docker build --no-cache -f 002-develop-yandasoft -t yandasoft .
docker tag yandasoft sord/yandasoft
docker push sord/yandasoft
