/* 
Copyright (c) 2023 Yuri Rubinsky (Chaosus)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
 */

#ifndef FLUIDSIM_CLASS_H
#define FLUIDSIM_CLASS_H

#ifdef WIN32
#include <windows.h>
#endif

#include <godot_cpp/classes/global_constants.hpp>
#include <godot_cpp/classes/image_texture.hpp>
#include <godot_cpp/classes/node.hpp>

#include <godot_cpp/core/binder_common.hpp>

using namespace godot;

class FluidSim : public Node {
	GDCLASS(FluidSim, Node);

	int N;
	double diffusion = 0.00001;
	double viscosity = 0.0000001;
	int iter = 4;

	float *s = nullptr;
	float *density = nullptr;

	float *Vx = nullptr;
	float *Vy = nullptr;

	float *Vx0 = nullptr;
	float *Vy0 = nullptr;

	bool initialized = false;

	Ref<Image> density_image;
	Ref<ImageTexture> density_texture;

	void set_bnd(int b, float *x);
	void lin_solve(int b, float *x, float *x0, float a, float c);
	void diffuse(int b, float *x, float *x0, float diff, float dt);
	void project(float *velocX, float *velocY, float *p, float *div);
	void advect(int b, float *d, float *d0, float *velocX, float *velocY, float dt);

protected:
	static void _bind_methods();

public:
	FluidSim();
	~FluidSim();

	void set_dimension(int p_size);
	int get_dimension() const;

	void set_diffusion(double p_diffusion);
	double get_diffusion() const;

	void set_viscosity(double p_viscosity);
	double get_viscosity() const;

	void set_iter(int p_iter);
	int get_iter() const;

	void add_density(const Vector2i &p_position, float p_amount);
	void add_velocity(const Vector2i &p_position, const Vector2 &p_amount);

	Ref<Texture> get_density_texture();

	void init();
	void step(float p_dt);
	void dispose();
};

#endif // FLUIDSIM_CLASS_H
