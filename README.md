# The fastest way to turn a program into a biobox

## How it works

### Container
A biobox is a (Docker) container which can contain a program and all its dependencies. bioboxes.org specifies a way how to pass and return data to such a container and standardizes common types of containers. A YAML file describes the container input and output. At first, you should pick a biobox type which you would like to implement.

### Writing tasks
Tasks are selected by the first string argument after `docker run biobox`. Simply place the command line which defines your task "mytask" in a text file named "mytask". POSIX shell syntax and all variables which you can list via `docker run bioboxes/base --list-var` are possible. For instance, additional arguments after the task definition are available via the environment variable `$BBX_ARGS`. Task names can also take the form of command line arguments, for instance `--list-var` is implemented as a task file.

### Adding programs
You should add native Debian packages from the official Debian repositories in the Dockerfile, if available. For other programs, create  folders like `opt/myprogram` in this folder. The path to these programs in the task definitions will be `$BBX_OPTDIR/myprogram`.

### Example
This template contains a working hello world biobox.

## What to do

1. Clone this GIT repository

        git clone https://github.com/bioboxes/biobox-template

2. Add required packages to `Dockerfile`. Put package names after `RUN dockerfile-install-packages`

3. Copy your software to `opt/`. You should create a subfolder under `opt/` in which you should place your scripts, programs and data.

4. Create a task file under `tasks/`. Name it like your task (spaces discouraged). Save it as a UTF8-encoded text file with UNIX line breaks (`\n`)

5. *Optional: edit `opt/default`*

6. *Optional: remove `tasks/helloworld` and `opt/helloworld_0.1/`*

7. Build the biobox

        docker build -t="dockeruser/name_of_biobox" .

8. Run the biobox

        docker run dockeruser/name_of_biobox name_of_task

9. *Optional: push the biobox to the Docker Hub. You need to register your docker user name (the prefix before /name_of_biobox) at hub.docker.com before uploading.*

        docker push dockeruser/name_of_biobox
