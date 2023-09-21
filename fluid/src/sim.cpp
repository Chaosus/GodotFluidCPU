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

#include "sim.h"

#include <godot_cpp/core/class_db.hpp>

#include <godot_cpp/classes/global_constants.hpp>
#include <godot_cpp/classes/label.hpp>
#include <godot_cpp/classes/multiplayer_api.hpp>
#include <godot_cpp/classes/multiplayer_peer.hpp>
#include <godot_cpp/variant/utility_functions.hpp>

using namespace godot;

#define IX(x, y) ((x) + (y)*N)

void FluidSim::_bind_methods() {
	ClassDB::bind_method(D_METHOD("set_dimension", "size"), &FluidSim::set_dimension);
	ClassDB::bind_method(D_METHOD("get_dimension"), &FluidSim::get_dimension);
	ClassDB::bind_method(D_METHOD("set_diffusion", "diffusion"), &FluidSim::set_diffusion);
	ClassDB::bind_method(D_METHOD("get_diffusion"), &FluidSim::get_diffusion);
	ClassDB::bind_method(D_METHOD("set_viscosity", "viscosity"), &FluidSim::set_viscosity);
	ClassDB::bind_method(D_METHOD("get_viscosity"), &FluidSim::get_viscosity);
	ClassDB::bind_method(D_METHOD("set_iter", "iter"), &FluidSim::set_iter);
	ClassDB::bind_method(D_METHOD("get_iter"), &FluidSim::get_iter);
	ClassDB::bind_method(D_METHOD("add_density", "position", "amount"), &FluidSim::add_density);
	ClassDB::bind_method(D_METHOD("add_velocity", "position", "amount"), &FluidSim::add_velocity);
	ClassDB::bind_method(D_METHOD("init"), &FluidSim::init);
	ClassDB::bind_method(D_METHOD("step", "delta"), &FluidSim::step);
	ClassDB::bind_method(D_METHOD("dispose"), &FluidSim::dispose);
	ClassDB::bind_method(D_METHOD("get_density_texture"), &FluidSim::get_density_texture);

	ADD_PROPERTY(PropertyInfo(Variant::INT, "dimension"), "set_dimension", "get_dimension");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "diffusion"), "set_diffusion", "get_diffusion");
	ADD_PROPERTY(PropertyInfo(Variant::FLOAT, "viscosity"), "set_viscosity", "get_viscosity");
	ADD_PROPERTY(PropertyInfo(Variant::INT, "iter"), "set_iter", "get_iter");
}

//////////////////////////////////////////////////////////////

void FluidSim::set_bnd(int b, float *x) {
	for (int i = 1; i < N - 1; i++) {
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
	}
	for (int j = 1; j < N - 1; j++) {
		x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
		x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
	}

	x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
	x[IX(0, N - 1)] = 0.5 * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
	x[IX(N - 1, 0)] = 0.5 * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
	x[IX(N - 1, N - 1)] = 0.5 * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
}

void FluidSim::lin_solve(int b, float *x, float *x0, float a, float c) {
	float cRecip = 1.0 / c;
	for (int t = 0; t < iter; t++) {
		for (int j = 1; j < N - 1; j++) {
			for (int i = 1; i < N - 1; i++) {
				x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip;
			}
		}
		set_bnd(b, x);
	}
}

void FluidSim::diffuse(int b, float *x, float *x0, float diff, float dt) {
	float a = dt * diff * (N - 2) * (N - 2);
	lin_solve(b, x, x0, a, 1 + 6 * a);
}

void FluidSim::project(float *velocX, float *velocY, float *p, float *div) {
	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			div[IX(i, j)] = -0.5f * (velocX[IX(i + 1, j)] - velocX[IX(i - 1, j)] + velocY[IX(i, j + 1)] - velocY[IX(i, j - 1)]) / N;
			p[IX(i, j)] = 0;
		}
	}

	set_bnd(0, div);
	set_bnd(0, p);
	lin_solve(0, p, div, 1, 6);

	for (int j = 1; j < N - 1; j++) {
		for (int i = 1; i < N - 1; i++) {
			velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) * N;
			velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) * N;
		}
	}

	set_bnd(1, velocX);
	set_bnd(2, velocY);
}

void FluidSim::advect(int b, float *d, float *d0, float *velocX, float *velocY, float dt) {
	float i0, i1, j0, j1;

	float dtx = dt * (N - 2);
	float dty = dt * (N - 2);

	float s0, s1, t0, t1;
	float tmp1, tmp2, x, y;

	float Nfloat = N - 2;
	float ifloat, jfloat;
	int i, j;

	for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
		for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
			tmp1 = dtx * velocX[IX(i, j)];
			tmp2 = dty * velocY[IX(i, j)];
			x = ifloat - tmp1;
			y = jfloat - tmp2;

			if (x < 0.5f)
				x = 0.5f;
			if (x > Nfloat + 0.5f)
				x = Nfloat + 0.5f;
			i0 = floorf(x);
			i1 = i0 + 1.0f;
			if (y < 0.5f)
				y = 0.5f;
			if (y > Nfloat + 0.5f)
				y = Nfloat + 0.5f;
			j0 = floorf(y);
			j1 = j0 + 1.0f;

			s1 = x - i0;
			s0 = 1.0f - s1;
			t1 = y - j0;
			t0 = 1.0f - t1;

			int i0i = i0;
			int i1i = i1;
			int j0i = j0;
			int j1i = j1;

			d[IX(i, j)] = s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) + s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
		}
	}
	set_bnd(b, d);
}

void FluidSim::step(float p_dt) {
	ERR_FAIL_COND(!initialized);

	diffuse(1, Vx0, Vx, viscosity, p_dt);
	diffuse(2, Vy0, Vy, viscosity, p_dt);

	project(Vx0, Vy0, Vx, Vy);

	advect(1, Vx, Vx0, Vx0, Vy0, p_dt);
	advect(2, Vy, Vy0, Vx0, Vy0, p_dt);

	project(Vx, Vy, Vx0, Vy0);
	diffuse(0, s, density, diffusion, p_dt);
	advect(0, density, s, Vx, Vy, p_dt);

	for (int y = 0; y < N; y++) {
		for (int x = 0; x < N; x++) {
			float val = density[IX(x, y)];
			density_image->set_pixel(x, y, Color(val, val, val));
		}
	}
	density_texture->update(density_image);
}

//////////////////////////////////////////////////////////////

void FluidSim::set_dimension(int p_size) {
	N = p_size;
}

int FluidSim::get_dimension() const {
	return N;
}

void FluidSim::set_diffusion(double p_diffusion) {
	diffusion = p_diffusion;
}

double FluidSim::get_diffusion() const {
	return diffusion;
}

void FluidSim::set_viscosity(double p_viscosity) {
	viscosity = p_viscosity;
}

double FluidSim::get_viscosity() const {
	return viscosity;
}

void FluidSim::set_iter(int p_iter) {
	iter = p_iter;
}

int FluidSim::get_iter() const {
	return iter;
}

void FluidSim::add_density(const Vector2i &p_position, float p_amount) {
	density[IX(p_position.x, p_position.y)] += p_amount;
}

void FluidSim::add_velocity(const Vector2i &p_position, const Vector2 &p_amount) {
	size_t index = IX(p_position.x, p_position.y);
	Vx[index] += p_amount.x;
	Vy[index] += p_amount.y;
}

Ref<Texture> FluidSim::get_density_texture() {
	return density_texture;
}

void FluidSim::init() {
	dispose();

	size_t sz = N * N * N;

	s = (float *)calloc(sz, sizeof(float));
	density = (float *)calloc(sz, sizeof(float));

	Vx = (float *)calloc(sz, sizeof(float));
	Vy = (float *)calloc(sz, sizeof(float));

	Vx0 = (float *)calloc(sz, sizeof(float));
	Vy0 = (float *)calloc(sz, sizeof(float));

	density_image.instantiate();
	density_image = Image::create(N, N, false, Image::FORMAT_RGB8);
	density_image->fill(Color(0, 0, 0));

	density_texture.instantiate();
	density_texture = ImageTexture::create_from_image(density_image);

	initialized = true;
}

void FluidSim::dispose() {
	if (initialized) {
		initialized = false;

		::free(s);
		::free(density);

		::free(Vx);
		::free(Vy);

		::free(Vx0);
		::free(Vy0);

		s = nullptr;
		density = nullptr;

		Vx = nullptr;
		Vy = nullptr;

		Vx0 = nullptr;
		Vy0 = nullptr;
	}
}

//////////////////////////////////////////////////////////////

FluidSim::FluidSim() {
}

FluidSim::~FluidSim() {
	dispose();
}