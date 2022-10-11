# Ray Tracer

Rust implementation for [The Ray Tracer Challenge][rt-challenge]. Download the CLI from
the [release][] page.

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

See comments in yaml scene files for information about how to get all resources (obj
models, ppm textures) used in the examples commands.

### Basic scene

The rendered scenes can be piped directly to [ImageMagick](https://imagemagick.org/) or
outputted to a file.

```
raytracer --scene samples/scenes/basic_scene.yaml | magick display
```

[![basic-scene](/samples/rendered/basic_scene.png?raw=true "basic-scene")](/samples/scenes/basic_scene.yaml)

### Books's cover

```
raytracer --scene samples/scenes/cover.yaml --out cover.ppm
```

[![cover](/samples/rendered/cover.png?raw=true "cover")](samples/scenes/cover.yaml)

### Spaceship model

```
raytracer --scene samples/scenes/space_ship.yaml \
          --obj samples/obj/space_ship.obj \
          --out space_ship.ppm
```

[![spaceship](/samples/rendered/space_ship.png?raw=true "spaceship")](samples/scenes/space_ship.yaml)

### Teapot and Spaceship models

Multiple `obj` models can be added/combined in a scene.

```
raytracer --scene samples/scenes/space_teapot.yaml \
          --obj samples/obj/space_ship.obj \
          --obj samples/obj/teapot_low.obj \
          --out space_teapot.ppm
```

[![space-teapot](/samples/rendered/space_teapot.png?raw=true "space-teapot")](samples/scenes/space_teapot.yaml)

### Constructive Solid Geometry (CSG)

```
raytracer --scene samples/scenes/csg.yaml --out csg.ppm
```

[![csg](/samples/rendered/csg.png?raw=true "csg")](samples/scenes/csg.yaml)

### Soft shadows (Area Light)

```
raytracer --scene samples/scenes/soft_shadows.yaml --out soft_shadows.ppm
```

[![soft-shadows](/samples/rendered/soft_shadows.png?raw=true "soft-shadows")](samples/scenes/soft_shadows.yaml)

### Bounding Boxes - Dragons

Beware that rendering this scene will be quite long (about 45 minutes on my machine with
16 cores).

```
raytracer --scene samples/scenes/dragons.yaml \
          --obj samples/obj/dragon.obj \
          --out dragons.ppm
```

[![dragons](/samples/rendered/dragons.png?raw=true "dragons")](samples/scenes/dragons.yaml)

### Texture Mapping - Patterns

```
raytracer --scene samples/scenes/checkered_cube.yaml --out checkered_cube.ppm
```

[![checkered-cube](/samples/rendered/checkered_cube.png?raw=true "checkered-cube")](samples/scenes/checkered_cube.yaml)

| [![checkered-sphere](/samples/rendered/checkered_sphere.png?raw=true "checkered-sphere")](samples/scenes/checkered_sphere.yaml) | [![checkered-plane](/samples/rendered/checkered_plane.png?raw=true "checkered-plane")](samples/scenes/checkered_plane.yaml) | [![checkered-cylinder](/samples/rendered/checkered_cylinder.png?raw=true "checkered-cylinder")](samples/scenes/checkered_cylinder.yaml) |
|---------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------|

### Texture Mapping - Image

```
raytracer --scene samples/scenes/earth.yaml \
          --ppm samples/textures/earthmap1k.ppm \
          --out earth.ppm
```

[![earth](/samples/rendered/earth.png?raw=true "earth")](samples/scenes/earth.yaml)

### Texture Mapping - Skybox

```
raytracer --scene samples/scenes/skybox.yaml \
          --ppm samples/textures/negx.ppm \
          --ppm samples/textures/negy.ppm \
          --ppm samples/textures/negz.ppm \
          --ppm samples/textures/posx.ppm \
          --ppm samples/textures/posy.ppm \
          --ppm samples/textures/posz.ppm \
          --out skybox.ppm
```

[![skybox](/samples/rendered/skybox.png?raw=true "skybox")](samples/scenes/skybox.yaml)
