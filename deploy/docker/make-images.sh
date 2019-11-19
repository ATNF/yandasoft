# Makes the base dependencies
docker build -f 001-yandabase -t yandabase .
# Makes yandasoft
docker build -f 002-develop-yandasoft -t yandasoft .
# Makes a slim version using the build products
docker build -f 003-yandaslim -t yandaslim .
