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

test:
	RUST_BACKTRACE=1 cargo test -- ${TESTFLAGS} ${ARGS}

$(TARGET): build

install: ${TARGET}
	sudo cp $? scripts/smolden /usr/bin/.
