# Ray Tracer

Rust implementation for [The Ray Tracer Challenge][rt-challenge]. Download the CLI from
the [release]() page.

[rt-challenge]: http://raytracerchallenge.com/
[release]: https://github.com/lerouxrgd/raytracer/releases

## Usage

```text
The Ray Tracer Challenge CLI

Usage: raytracer [OPTIONS] --scene <FILE>

Options:
      --scene <FILE>  A yaml description of the scene to render
      --obj <FILE>    Optional obj models to add to the scene
      --ppm <FILE>    Optional ppm textures to use as material
      --out <OUT>     Optional output ppm file, defaults to stdout
  -h, --help          Print help information
```

## Examples

Check out comments in yaml scene files to get all resources (obj models, ppm textures)
used in the examples commands.

The rendered scene can be piped directly to [ImageMagick](https://imagemagick.org/) or
outputted to a file.

### Basic scene

```
raytracer --scene samples/scenes/basic_scene.yaml | magick display
```

![basic-scene](/samples/rendered/basic_scene.png?raw=true "basic-scene")

### Books's cover

```
raytracer --scene samples/scenes/cover.yaml --out cover.ppm
```

![cover](/samples/rendered/cover.png?raw=true "cover")

### Spaceship model

```
raytracer --scene samples/scenes/space_ship.yaml \
          --obj samples/obj/space_ship.obj \
          --out space_ship.ppm
```

![spaceship](/samples/rendered/space_ship.png?raw=true "spaceship")

### Teapot and Spaceship models

Multiple `obj` models can be added/combined in a scene.

```
raytracer --scene samples/scenes/space_teapot.yaml \
          --obj samples/obj/space_ship.obj \
          --obj samples/obj/teapot_low.obj \
          --out space_teapot.ppm
```

![space-teapot](/samples/rendered/space_teapot.png?raw=true "space-teapot")
