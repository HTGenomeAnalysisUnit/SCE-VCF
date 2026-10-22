FROM debian:bookworm-slim AS builder

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        git \
        libhts-dev \
        nim \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY sceVCF.nimble ./
COPY src/ ./src/

RUN nimble install -y --depsOnly \
    && nimble build -y -d:release --opt:size --passL:-s \
    && install -Dm755 bin/sceVCF /out/sceVCF

FROM debian:bookworm-slim AS runtime

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        libhts3 \
        libpcre3 \
    && rm -rf /var/lib/apt/lists/* \
    # hts-nim loads libhts.so, which is normally supplied by libhts-dev.
    && ln -s "$(ldconfig -p | awk '$1 == "libhts.so.3" { print $NF; exit }')" /usr/local/lib/libhts.so \
    && ldconfig

COPY --from=builder /out/sceVCF /usr/local/bin/sceVCF
COPY LICENSE /usr/share/doc/sceVCF/LICENSE

# Check the tool runs in the final image
RUN sceVCF --help

WORKDIR /data
ENTRYPOINT ["/usr/local/bin/sceVCF"]
CMD ["--help"]
