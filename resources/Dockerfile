FROM ghcr.io/prefix-dev/pixi:bookworm-slim

WORKDIR /app

COPY ./pixi.toml .

RUN pixi global install && rm -rf ~/.cache/rattler