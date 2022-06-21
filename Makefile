BASE = /home/brent/Projects/scope
ARGS =
TESTFLAGS = --test-threads=1 --nocapture
SHORT = 0

TARGET=/home/brent/Projects/scope/target/x86_64-unknown-linux-gnu/release/scope

ifeq ($(SHORT),0)
TESTFLAGS += --include-ignored
endif

WOODS_DEST = 'woods:Programs/scope/scope-bin'

build: src/*.rs
	RUSTFLAGS='-C target-feature=+crt-static' \
		cargo build --release --target x86_64-unknown-linux-gnu \

woods: build
	scp -C ${TARGET} ${WOODS_DEST}

# ELAND_DEST = 'eland:programs/semp/.'
# eland:
# # see https://msfjarvis.dev/posts/building-static-rust-binaries-for-linux
# 	RUSTFLAGS='-C target-feature=+crt-static' \
# 	cargo build --release --target x86_64-unknown-linux-gnu
# 	scp -C ${BASE}/target/x86_64-unknown-linux-gnu/release/rust-semp ${ELAND_DEST}

# eland.scripts: scripts/time.awk
# 	scp -C $? ${ELAND_DEST}

test:
	RUST_BACKTRACE=1 cargo test -- ${TESTFLAGS} ${ARGS}

install: ${TARGET}
	cp $? /usr/bin/.

# profile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
# 	valgrind --tool=callgrind --callgrind-out-file=callgrind.out	\
# 		--collect-jumps=yes --simulate-cache=yes		\
# 		${BASE}/target/release/$(1)

# profile.one_iter:
# 	$(call profile,one_iter)

# profile.num_jac:
# 	$(call profile,num_jac)

# profile.full:
# 	$(call profile,full)

# memprofile = RUSTFLAGS='-g' cargo build --release --bin $(1); \
# 		heaptrack ${BASE}/target/release/$(1)

# memprofile.full:
# 	$(call memprofile,full)

# time:
# 	RUSTFLAGS='-g' cargo build --release --bin full
# 	sh -c "time ${BASE}/target/release/full"
