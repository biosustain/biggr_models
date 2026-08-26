bigg_models
-----------


[BiGG](https://bigg.bio) is a website for browsing normalized, gold-standard genome-scale models and the corresponding metabolite and reaction namespace.
BiGG is based on [legacy BiGG](http://bigg.ucsd.edu), described in the following publication:

King ZA, Lu JS, Dräger A, Miller PC, Federowicz S, Lerman JA, Ebrahim A, Palsson BO, and Lewis NE. (2015). BiGG Models: A platform for integrating, standardizing, and sharing genome-scale models. Nucl Acids Res. doi:[10.1093/nar/gkv1049](https://doi.org/10.1093/nar/gkv1049).

This repository includes the web server and front-end for BiGG. The database is managed by [COBRAdb](https://github.com/biosustain/cobradb), Escher maps are auto-generated with the help of [BiGG Maps](http://github.com/biosustain/bigg_maps), and the [BiGG Python library](http://github.com/biosustain/bigg) can be used to easily access the database and its API.

Installation
============

It is recommended to install and run BiGG using [Docker](https://www.docker.com/).

To run the server using Docker, make sure Docker is installed and then follow these steps:

For development environments:
1. Download cobradb with ```git clone git@github.com:biosustain/cobradb.git```.
2. Download the code with ```git clone git@github.com:biosustain/bigg_models.git```
3. ```cd bigg_models```
4. For setting up the `settings.ini` and `.env` files, and various `/data/...` directories,
    please refer to the [cobradb](https://github.com/biosustain/cobradb) readme.
5. Run ```docker compose --profile dev up --build```.

For production environments:
1. Download the code with ```git clone git@github.com:biosustain/bigg_models.git```
2. ```cd bigg_models```
3. For setting up the `settings.ini` and `.env` files, and various `/data/...` directories,
    please refer to the [cobradb](https://github.com/biosustain/cobradb) readme.
4. Run ```docker compose --profile prod up --build```.

Note that for the bigg.bio server, a CI/CD pipeline is implemented on GitHub.

Redeployment
============

The production stack runs on this [VM](https://portal.azure.com/#@dtudk.onmicrosoft.com/resource/subscriptions/aee8556f-d2fd-4efd-a6bd-f341a90fa76e/resourceGroups/RG-RECON/providers/Microsoft.Compute/virtualMachines/BIGGR-PROD/overview) at ```~/projects/biggr_models```.

**Deploying new code (recommended).** Push to `main` on GitHub. The
[deploy-prod](.github/workflows/deploy-prod.yml) workflow copies the code to the
VM and recreates the containers. This is the preferred route for any change to
the application code.

**Restarting or redeploying manually.** If the containers are down, or if you
need to redeploy without a code change, SSH into the VM and run:

    cd /data/projects/biggr_models/
    docker compose --profile prod down
    docker compose --profile prod up -d --build --force-recreate --remove-orphans

Check the result with
`docker compose ps` and `docker compose logs -f biggr-web`.

**About automatic restarts.** The `biggr-web`, `nginx-prod`, and `biggr-db`
services are configured with `restart: always`, so they come back up after a
crash or a VM reboot. Note that this does not fix the underlying cause: if the
disk is full or a newly deployed build is broken, the container will simply
restart in a loop. In those cases you still need to inspect the VM and push a fix.

**Renewing SSL certificates.** To renew the certificates manually, run:

    docker compose --profile prod up certbot && docker compose --profile prod restart nginx-prod
    