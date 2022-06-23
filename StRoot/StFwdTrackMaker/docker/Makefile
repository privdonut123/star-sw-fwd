
NODEBUG=yes
STAR=/star-sw
STAR_LIB=$STAR/.${STAR_HOST_SYS}/LIB
STAR_BIN=$STAR/.${STAR_HOST_SYS}/BIN
STAR_SYS=x8664_sl7
PATH=$STAR_BIN:$STAR/mgr:$PATH
LD_LIBRARY_PATH=$STAR_LIB:$LD_LIBRARY_PATH
ROOT_INCLUDE_PATH=$STAR/.${STAR_HOST_SYS}/include



# HOST start docker container
start:
	docker run --rm -ti --network none -v /Users/jdb/bnl/work/ssw/tracking:/work/ --name fwd-tracking star-spack:root5

# HOST attach to docker container
a: attach

# HOST attach to docker container
attach:
	docker attach fwd-tracking
