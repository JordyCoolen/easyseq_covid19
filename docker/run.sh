# ${1} <runtimecontainername> example: covid
# ${2} <imagename> example: jonovox/easyseq_covid19:latest

docker run \
  -it \
  --rm \
  --name ${1} \
  --mount type=bind,source=${PWD},target=/workflow \
  ${2} /bin/bash