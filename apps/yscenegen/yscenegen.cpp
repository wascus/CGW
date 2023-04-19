//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include <yocto/yocto_color.h>
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>

using namespace yocto;

#include <filesystem>
#include <memory>

#include "ext/perlin-noise/noise1234.h"

float noise(const vec3f& p) { return noise3(p.x, p.y, p.z); }
vec2f noise2(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11})};
}
vec3f noise3(const vec3f& p) {
  return {noise(p + vec3f{0, 0, 0}), noise(p + vec3f{3, 7, 11}),
      noise(p + vec3f{13, 17, 19})};
}
float fbm(const vec3f& p, int octaves) {
    rand();
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float turbulence(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 1.0f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * fabs(noise(p * scale));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}
float ridge(const vec3f& p, int octaves) {
  auto sum    = 0.0f;
  auto weight = 0.5f;
  auto scale  = 1.0f;
  for (auto octave = 0; octave < octaves; octave++) {
    sum += weight * (1 - fabs(noise(p * scale))) * (1 - fabs(noise(p * scale)));
    weight /= 2;
    scale *= 2;
  }
  return sum;
}

sceneio_instance* get_instance(sceneio_scene* scene, const string& name) {
  for (auto instance : scene->instances)
    if (instance->name == name) return instance;
  print_fatal("unknown instance " + name);
  return nullptr;
}

void add_polyline(sceneio_shape* shape, const vector<vec3f>& positions,
    const vector<vec4f>& colors, float thickness = 0.0001f) {
  auto offset = (int)shape->positions.size();
  shape->positions.insert(
      shape->positions.end(), positions.begin(), positions.end());
  shape->colors.insert(shape->colors.end(), colors.begin(), colors.end());
  shape->radius.insert(shape->radius.end(), positions.size(), thickness);
  for (auto idx = 0; idx < positions.size() - 1; idx++) {
    shape->lines.push_back({offset + idx, offset + idx + 1});
  }
}

void sample_shape(vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, sceneio_shape* shape, int num) {
  auto triangles  = shape->triangles;
  auto qtriangles = quads_to_triangles(shape->quads);
  triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
  auto cdf = sample_triangles_cdf(triangles, shape->positions);
  auto rng = make_rng(19873991);
  for (auto idx = 0; idx < num; idx++) {
    auto [elem, uv] = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    auto q          = triangles[elem];
    positions.push_back(interpolate_triangle(shape->positions[q.x],
        shape->positions[q.y], shape->positions[q.z], uv));
    normals.push_back(normalize(interpolate_triangle(
        shape->normals[q.x], shape->normals[q.y], shape->normals[q.z], uv)));
    if (!texcoords.empty()) {
      texcoords.push_back(interpolate_triangle(shape->texcoords[q.x],
          shape->texcoords[q.y], shape->texcoords[q.z], uv));
    } else {
      texcoords.push_back(uv);
    }
  }
}

struct terrain_params {
  float size    = 0.1f;
  vec3f center  = zero3f;
  float height  = 0.1f;
  float scale   = 10;
  int   octaves = 8;
  vec4f bottom  = srgb_to_rgb(vec4f{154, 205, 50, 255} / 255);
  vec4f middle  = srgb_to_rgb(vec4f{205, 133, 63, 255} / 255);
  vec4f top     = srgb_to_rgb(vec4f{240, 255, 255, 255} / 255);
};




auto rng = make_rng(71);

float distance3D(const vec3f p1, const vec3f p2) {
  return sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2) + pow(p2.z - p1.z, 2));
}

void make_terrain(sceneio_scene* scene, sceneio_instance* instance,
    const terrain_params& params) {
  // YOUR CODE GOES HERE ---------------------------------------------------
    float max_y = 0;
  for (auto& n : instance->shape->positions) {
    auto a = (1 - distance3D(n, params.center) / params.size);
    auto r = ridge(n * params.scale, params.octaves);
    n.y    = (r * a) * params.size;
    max_y = (n.y > max_y) ? n.y : max_y;

  }
    for (auto i = 0; i < instance->shape->positions.size(); i++)
    {
        auto n = instance->shape->positions[i];
        if (n.y > 0.60  * params.height)
            instance->shape->colors.push_back(params.top);
        else if (n.y > 0.30 * params.height)
            instance->shape->colors.push_back(params.middle);
        else
            instance->shape->colors.push_back(params.bottom);
    }
  instance->shape->normals = compute_normals(
      instance->shape->quads, instance->shape->positions);
}

struct displacement_params {
  float height  = 0.02f;
  float scale   = 50;
  int   octaves = 8;
  vec4f bottom  = srgb_to_rgb(vec4f{64, 224, 208, 255} / 255);
  vec4f top     = srgb_to_rgb(vec4f{244, 164, 96, 255} / 255);
};

vec3f translate(vec3f p, vec3f v, float t) {
  auto n = normalize(v) * t;
  return p + n;
}

void make_displacement(sceneio_scene* scene, sceneio_instance* instance,
    const displacement_params& params) {
  // YOUR CODE GOES HERE ---------------------------------------------------
  int   c          = 0;
  float max_height = 0;
  for (auto& pos : instance->shape->positions) {
    auto t = turbulence(pos * params.scale, params.octaves) * params.height;
    auto pos_copy = pos;
    pos           = translate(pos, instance->shape->normals[c++], t);
    if (distance3D(pos, pos_copy) > max_height)
      max_height = distance3D(pos, pos_copy);

    instance->shape->colors.push_back((interpolate_line(
        params.bottom, params.top, distance3D(pos, pos_copy) / params.height)));
  }
  instance->shape->normals = compute_normals(
      instance->shape->quads, instance->shape->positions);
}

struct hair_params {
  int   num      = 10000;
  int   steps    = 1;
  float lenght   = 0.02f;
  float scale    = 250;
  float strength = 0.01f;
  float gravity  = 0.0f;
  vec4f bottom   = srgb_to_rgb(vec4f{25, 25, 25, 255} / 255);
  vec4f top      = srgb_to_rgb(vec4f{244, 164, 96, 255} / 255);
};

void make_hair(sceneio_scene* scene, sceneio_instance* instance,
    sceneio_instance* hair, const hair_params& params) {
    hair->shape = add_shape(scene, "hair");
    auto num = params.num*10;
    // YOUR CODE GOES HERE ---------------------------------------------------
    int initial_position_size = (int)instance->shape->positions.size();
    sample_shape(instance->shape->positions, instance->shape->normals,
        instance->shape->texcoords, instance->shape, num);



    int i = 0;

    for (i = initial_position_size; i < num + initial_position_size ; i++)
    {
        vector<vec3f> positions = {};
        positions.push_back(instance->shape->positions[i]);
        vector<vec4f> colors = {};
        vec3f normal = instance->shape->normals[i];
        vec3f pos = instance->shape->positions[i];;
        vec3f pos2 = zero3f;
        
        for (auto j = 0; j < params.steps; j++)
        {
            pos2 = pos + normalize(normal)*(params.lenght/params.steps) + noise3(pos*params.scale)*params.strength;
            pos2.y -= params.gravity;
            positions.push_back(pos2);
            auto c = interpolate_line(params.bottom, params.top, (float)j/params.steps);
            colors.push_back(c);
            normal = pos2 - pos;
            pos = pos2;
        }
        colors.push_back(params.top);

        add_polyline(hair->shape, positions, colors);

    }
}

struct grass_params {
  int num = 10000;
};

int randN(int n)
/* Genera interi casuali da 0 a n*/
{
    return (int)(rand1f(rng)*n);
}

int randInt(int a, int b, rng_state rng)
/* Genera interi casuali in [a, b]*/
{
    return (int)(rand1f(rng)*(b - a)) + a;
}

float randFloat(float a, float b, rng_state rng)
{
    return (rand1f(rng)*(b - a)) + a;
}

vec3f v3f(float f)
{
    return vec3f({f, f, f});
}

namespace yocto {
 frame3f my_rotation_frame(const vec3f& axis, float angle) {
  auto s = sin(angle), c = cos(angle);
  auto vv = normalize(axis);
  return {{c + (axis.x - c) * vv.x * vv.x, (axis.x - c) * vv.x * vv.y + s * vv.z,
              (axis.x - c) * vv.x * vv.z - s * vv.y},
      {(axis.y - c) * vv.x * vv.y - s * vv.z, c + (axis.y - c) * vv.y * vv.y,
          (axis.y - c) * vv.y * vv.z + s * vv.x},
      {(axis.z - c) * vv.x * vv.z + s * vv.y, (axis.z - c) * vv.y * vv.z - s * vv.x,
          c + (axis.z - c) * vv.z * vv.z},
      {0, 0, 0}};
 }
}

void make_grass(sceneio_scene* scene, sceneio_instance* object,
    const vector<sceneio_instance*>& grasses, const grass_params& params) {
    // YOUR CODE GOES HERE ---------------------------------------------------
    auto randy = make_rng(0);
    auto initial_position_size = (int)object->shape->positions.size();
    sample_shape(object->shape->positions, object->shape->normals, object->shape->texcoords, object->shape, params.num);
    auto i = 0;
    string name = "grass";

    for ( i = initial_position_size; i < initial_position_size + params.num; i++)
    {
        sceneio_instance *nn = add_instance(scene, name.append(std::to_string(i)));
        int r = randN((int)grasses.size());
        nn->shape = grasses[r]->shape;
        nn->material = grasses[r]->material;


        nn->frame = nn->frame * translation_frame(object->shape->positions[i]);

        auto rotY = rand1f(randy)*2*pi;
        auto rotZ = rand1f(randy)*0.1 + 0.1;

        auto s = vec3f{rand1f(randy)/10.0f + 0.9f, rand1f(randy)/10.0f + 0.9f, rand1f(randy)/10.0f + 0.9f};

        nn->frame *= scaling_frame(s) * rotation_frame({0, 1, 0}, rotY) * rotation_frame({0, 0, 1}, rotZ);

  
        
    }

}



bool between(float m, float a, float b)
{
    return m >= a && m <= b;
}

vector<vec3f> generatePoints(float xmin, float xmax, float ymin, float ymax, int n, float k, rng_state rng)
{
    vector<vec3f> res = {};
    float square_length = interpolate_line(ymin, ymax, (float)1/(float)n) - ymin;
    float s = square_length/1.9f;
    for (int i = 1; i <= n; i++)
    {
        float y = interpolate_line(ymin, ymax, (float)i/(float)n) - square_length/(float)2.0;
        for (int j = 1; j <= n; j++)
        {
            float x = interpolate_line(xmin, xmax, (float)j/(float)n) - square_length/(float)2.0;
            res.push_back(vec3f{x + (rand1f(rng)*s - s/2.0f), 0, y + (rand1f(rng)*s - s/2.0f)});
        }
    }
    return res;
}

struct extra1_params   
{
    vec4f plane_color = vec4f{1, 1, 1, 1};
    vec4f points_color = vec4f{0, 0, 0, 1};
    float k = 0;
    int num = 20;
    

};

void make_uniformed_points(sceneio_scene* scene, sceneio_instance* instance, const extra1_params& params)
{
    int positions_size = (int)instance->shape->positions.size();
    auto randy = make_rng(0);
    auto colors = instance->shape->colors;
    for (int i = 0; i < positions_size; i++)
    {
        instance->shape->colors.push_back(params.plane_color);
    }

    float xmin = 0;
    float xmax = 0;
    float zmin = 0;
    float zmax = 0;

    for (int i = 0; i < instance->shape->positions.size(); i++)
    {
        auto pos = instance->shape->positions[i];
        if (pos.x < xmin) xmin = pos.x;
        if (pos.x > xmax) xmax = pos.x;
        if (pos.z < zmin) zmin = pos.z;
        if (pos.z > zmax) zmax = pos.z;
    }

    positions_size = (int)instance->shape->positions.size();
    vec3f center = {(xmin + xmax)/2, 0, (zmin + zmax)/2};
    float width = xmax - xmin;
    float radius = width/100;

    vector<vec3f> points = generatePoints(xmin, xmax, zmin, zmax, params.num, params.k, randy);

    for (auto& point : points)
    {
        for  (int i = 0; i < positions_size; i++)
        {
            if (distance(instance->shape->positions[i], point) <= radius/2.0f)
                instance->shape->colors[i] = params.points_color;
        }
    }
}

struct extra2_params
{
    int num = 10'000;
    int small = 3500;
    int medium = 2000;
    int big = 100;
    vec4f bottom =  srgb_to_rgb(vec4f{0, 0, 0, 255} / 255);  
    vec4f top =   srgb_to_rgb(vec4f{82, 82, 82, 255} / 255);  
};

void make_lunar_surface(sceneio_scene* scene, sceneio_instance* instance,
    const extra2_params& params)
{
    auto dp = displacement_params{  0.02f, 1, 8, srgb_to_rgb(vec4f{120,120,120, 255} / 255), srgb_to_rgb(vec4f{220,220,220, 255} / 255)};
    make_displacement(scene, instance, dp);
    auto randy = make_rng(1);
    int initial_size = (int)instance->shape->positions.size();
    auto positions_copy = instance->shape->positions;
    sample_shape(instance->shape->positions, instance->shape->normals, instance->shape->texcoords, instance->shape, params.num);

    // adding big holes
    for (int ii = 0; ii < params.big; ii++)
    {
        int i = randInt(initial_size, initial_size + params.num, randy);
        auto radius = rand1f(randy)*0.01;
        for (int j = 0; j < (int)positions_copy.size(); j++)
            if (distance3D(positions_copy[j], instance->shape->positions[i]) <= radius)
            {
                if (distance3D(positions_copy[j], instance->shape->positions[i]) > radius - 0.001)
                    instance->shape->positions[j] =  positions_copy[j] + normalize(instance->shape->normals[j])*0.0005;
                else
                    instance->shape->positions[j] =  positions_copy[j] - normalize(instance->shape->normals[j])*0.001;
            }
    }

    // adding medium holes
    for (int ii = 0; ii < params.medium; ii++)
    {
        int i = randInt(initial_size, initial_size + params.num, randy);
        auto radius = rand1f(randy)*0.008;
        for (int j = 0; j < (int)positions_copy.size(); j++)
            if (distance3D(positions_copy[j], instance->shape->positions[i]) <= radius) 
            {
                if (distance3D(positions_copy[j], instance->shape->positions[i]) > radius - 0.001)
                    instance->shape->positions[j] =  positions_copy[j] + normalize(instance->shape->normals[j])*0.0005;
                else
                    instance->shape->positions[j] =  positions_copy[j] - normalize(instance->shape->normals[j])*0.001;
            }
    }


    // adding small holes
    for (int ii = 0; ii < params.small; ii++)
    {
        auto i = randInt(initial_size, initial_size + params.num, randy);
        auto radius = rand1f(randy)*0.0008;
        for (int j = 0; j < (int)positions_copy.size(); j++)
            if (distance3D(positions_copy[j], instance->shape->positions[i]) <= radius) 
                instance->shape->positions[j] =  positions_copy[j] - normalize(instance->shape->normals[j])*0.003;

    }

    for (int i = 0; i < instance->shape->positions.size(); i++)
        instance->shape->colors.push_back(interpolate_line(params.bottom, params.top, distance3D(instance->shape->positions[i], zero3f)/0.02f));
    instance->shape->normals = compute_normals(
        instance->shape->quads, instance->shape->positions);


}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto terrain      = ""s;
  auto tparams      = terrain_params{};
  auto displacement = ""s;
  auto dparams      = displacement_params{};
  auto hair         = ""s;
  auto hairbase     = ""s;
  auto hparams      = hair_params{};
  auto grass        = ""s;
  auto grassbase    = ""s;
  auto gparams      = grass_params{};
  auto output       = "out.json"s;
  auto filename     = "scene.json"s;

  auto extra1       = ""s;
  auto extra2       = ""s;
  auto e1params = extra1_params{};
  auto e2params = extra2_params{};
  // parse command line
  auto cli = make_cli("yscenegen", "Make procedural scenes");
  add_option(cli, "--terrain", terrain, "terrain object");
  add_option(cli, "--displacement", displacement, "displacement object");
  add_option(cli, "--hair", hair, "hair object");
  add_option(cli, "--hairbase", hairbase, "hairbase object");
  add_option(cli, "--grass", grass, "grass object");
  add_option(cli, "--grassbase", grassbase, "grassbase object");
  add_option(cli, "--hairnum", hparams.num, "hair number");
  add_option(cli, "--hairlen", hparams.lenght, "hair length");
  add_option(cli, "--hairstr", hparams.strength, "hair strength");
  add_option(cli, "--hairgrav", hparams.gravity, "hair gravity");
  add_option(cli, "--hairstep", hparams.steps, "hair steps");
  add_option(cli, "--output,-o", output, "output scene");
  add_option(cli, "scene", filename, "input scene", true);
  ///////////
  add_option(cli, "--extra1", extra1, "extra object");
  add_option(cli, "--extra2", extra2, "extra2 object");

  parse_cli(cli, argc, argv);

  // load scene
  auto scene_guard = std::make_unique<sceneio_scene>();
  auto scene       = scene_guard.get();
  auto ioerror     = ""s;
  if (!load_scene(filename, scene, ioerror, print_progress))
    print_fatal(ioerror);

  // create procedural geometry
  if (terrain != "") {
    make_terrain(scene, get_instance(scene, terrain), tparams);
  }
  if (displacement != "") {
    make_displacement(scene, get_instance(scene, displacement), dparams);
  }
  if (hair != "") {
    make_hair(scene, get_instance(scene, hairbase), get_instance(scene, hair),
        hparams);
  }
  if (grass != "") {
    auto grasses = vector<sceneio_instance*>{};
    for (auto instance : scene->instances)
      if (instance->name.find(grass) != scene->name.npos)
        grasses.push_back(instance);
    make_grass(scene, get_instance(scene, grassbase), grasses, gparams);
  }

  if (extra1 != "")
  {
      make_uniformed_points(scene, get_instance(scene, extra1), e1params);
  }
  if (extra2 != "")
  {
    make_lunar_surface(scene, get_instance(scene, extra2), e2params);
  }

  // make a directory if needed
  if (!make_directory(path_dirname(output), ioerror)) print_fatal(ioerror);
  if (!scene->shapes.empty()) {
    if (!make_directory(path_join(path_dirname(output), "shapes"), ioerror))
      print_fatal(ioerror);
  }
  if (!scene->textures.empty()) {
    if (!make_directory(path_join(path_dirname(output), "textures"), ioerror))
      print_fatal(ioerror);
  }

  // save scene
  if (!save_scene(output, scene, ioerror, print_progress)) print_fatal(ioerror);

  // done
  return 0;
}
