#Compile a static build with no dependencies on hts lib
#See: https://github.com/brentp/hts-nim#static-builds
#docker://brentp/musl-hts-nim:latest

pkg_name="sceVCF"

singularity exec \
--bind $PWD \
--bind $PWD:/load/ \
--bind /project/alfredo \
/project/alfredo/singularity/musl-hts-nim_latest.sif \
/usr/local/bin/nsb \
-s $PWD/src/${pkg_name}.nim \
--nimble-file $PWD/${pkg_name}.nimble -- -d:danger -d:release
