#! /usr/bin/env python

import shutil
import subprocess
import os
from harpy.common.conda import CONDA_ENVS

dockerfile_text = """
FROM ghcr.io/prefix-dev/pixi:0.63.2 AS build

# copy source code, pixi.toml and pixi.lock to the container
WORKDIR /app
COPY . .

# use `--locked` to ensure the lockfile is up to date with pixi.toml
RUN pixi install --locked && rm -rf ~/.cache/rattler

# create the shell-hook bash script to activate the environment
RUN echo "#!/bin/bash" > /app/entrypoint.sh && \\
    pixi shell-hook -s bash >> /app/entrypoint.sh && \\
    echo 'exec "$@"' >> /app/entrypoint.sh && \\
    chmod +x /app/entrypoint.sh

FROM ubuntu:24.04 AS production
WORKDIR /app
COPY --from=build --chmod=0755 /app/entrypoint.sh /app/entrypoint.sh
COPY --from=build /app/.pixi/envs/default /app/.pixi/envs/default

ENTRYPOINT ["/app/entrypoint.sh"]
"""

def create_pixi_dockerfiles():
    '''
    Using the defined environments, create a series of folders where each has a dockerfile
    and pixi.toml file to create one of the environments.
    '''
    shutil.rmtree("container", ignore_errors=True)
    for env,deps in CONDA_ENVS.items():
        os.makedirs(f"container/{env}", exist_ok=True)
        with open(f"container/{env}/Dockerfile", "w") as dockerfile:
            dockerfile.write(dockerfile_text)
        if env == "report":
            subprocess.run(
                f"pixi init container/{env} -c conda-forge -c r".split(),
                check = True
            )
        else:
            subprocess.run(
                f"pixi init container/{env} -c conda-forge -c bioconda".split(),
                check = True    
            )

        subprocess.run(
            ["pixi", "add", "--no-progress", "--manifest-path", f"container/{env}/pixi.toml", *deps],
            check = True
        )
        shutil.rmtree("container/.pixi", ignore_errors=True)
