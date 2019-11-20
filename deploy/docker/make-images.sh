# Makes the base dependencies
docker build --no-cache -f 001-yandabase -t yandabase .
# Makes yandasoft
docker build --no-cache -f 002-develop-yandasoft -t yandasoft .
# Makes a slim version using the build products
docker build --no-cache -f 003-yandaslim -t yandaslim .
