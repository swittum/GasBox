#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <random>
#include <iomanip>
#include "hdf5.h"

using std::cout;
using std::endl;

typedef std::vector<double> vd;
typedef std::vector<int> vi;
typedef std::mt19937 rn;

double skp(vd v1, vd v2)
{
    double out = 0.;
    for (int i = 0; i < v1.size(); i++)
    {
        out += v1[i] * v2[i];
    }
    return out;
}

double norm(vd v1)
{
    double norm2 = skp(v1, v1);
    double out = sqrt(norm2);
    return out;
}

vd substract(vd v1, vd v2)
{
    vd out(v1.size());
    for (int i = 0; i < v1.size(); i++)
    {
        out[i] = v1[i] - v2[i];
    }
    return out;
}

vd scale(vd v, double lmd)
{
    vd out(v.size());
    for (int i = 0; i < v.size(); i++)
    {
        out[i] = lmd*v[i];
    }
    return out;
}

class Particle
{
public:
    double mass, radius;
    vd x, v;

    Particle(double m, double r)
    {
        mass = m;
        radius = r;
        init();
    }

    ~Particle()
    {
        // cout << "Destroying member of Particle class" << endl;
    }

    void move(double dt)
    {
        for (int i = 0; i < 2; i++)
        {
            x[i] += dt * v[i];
        }
    }

    double get_v_abs()
    {
        return norm(v);
    }

private:
    void init()
    {
        std::random_device rd;
        rn gen(rd());
        for (int i = 0; i < 2; i++)
        {
            x.push_back(gen() / static_cast<double>(gen.max()) * (1 - 2 * radius) + radius);
            v.push_back(gen() / static_cast<double>(gen.max()) - 0.5);
        }
    }
};

class GasContainer
{
public:
    std::vector<Particle> particles;
    int n_particles;
    double dt, radius, mass;
    std::vector<vi> combinations;
    vd distances;
    std::vector<std::vector<vd> > pos_animation;
    std::vector<vd> vel_animation;

    GasContainer(int n, double m, double r, double h)
    {
        n_particles = n;
        dt = h;
        radius = r;
        mass = m;
        compute_combinations();
        for (int i = 0; i < n_particles; i++)
        {
            Particle p(mass, radius);
            particles.push_back(p);
        }
    }

    void evolve(double t_end)
    {
        int index_end = round(t_end / dt);
        std::vector<vd> circles;
        vd vels;

        for (int k = 0; k < index_end; k++)
        {
            for (int i = 0; i < n_particles; i++)
            {
                particles[i].move(dt);
            }
            wall_collision();
            particle_collision();
            vels.clear();
            circles.clear();
            for (int i = 0; i < n_particles; i++)
            {
                circles.push_back(particles[i].x);
                vels.push_back(particles[i].get_v_abs());
            }
            vel_animation.push_back(vels);
            pos_animation.push_back(circles);
        }
        cout << "Running" << endl;
    }

    void wall_collision()
    {
        vd pos, vel;
        for (int i = 0; i < n_particles; i++)
        {
            pos = particles[i].x;
            for (int j = 0; j < 2; j++)
            {
                if (pos[j] < radius || pos[j] > 1 - radius)
                {
                    particles[i].v[j] *= -1;
                }
            }
        }
    }

    void particle_collision()
    {
        int n = combinations.size();
        vd r1, r2, delta_r, v_tmp1, v_tmp2;
        vi inds;
        for (int i = 0; i < n; i++)
        {
            inds = combinations[i];
            r1 = particles[inds[0]].x;
            r2 = particles[inds[1]].x;
            delta_r.clear();
            for (int k = 0; k < r1.size(); k++)
            {
                delta_r.push_back(r1[k]-r2[k]);
            }
            if (norm(delta_r) < 2 * radius)
            {
                v_tmp1 = scattering(particles[combinations[i][0]], particles[combinations[i][1]]);
                v_tmp2 = scattering(particles[combinations[i][1]], particles[combinations[i][0]]);
                particles[combinations[i][0]].v = v_tmp1;
                particles[combinations[i][1]].v = v_tmp2;
            }
        }
    }

    vd scattering(Particle p1, Particle p2)
    {
        vd v_new;
        double coef = skp(substract(p1.v, p2.v), substract(p1.x, p2.x)) / pow(norm(substract(p1.x, p2.x)), 2);
        for (int i = 0; i < 2; i++)
        {
            v_new.push_back(p1.v[i] - scale(substract(p1.x, p2.x), coef)[i]);
        }
        return v_new;
    }

    void compute_combinations()
    {
        vi combi;
        for (int i = 0; i < n_particles; i++)
        {
            for (int j = i + 1; j < n_particles; j++)
            {
                combi.clear();
                combi.push_back(i);
                combi.push_back(j);
                combinations.push_back(combi);
            }
        }
    }

    ~GasContainer()
    {
        cout << "Destroying member of GasContainer class" << endl;
    }
};

int main()
{
    cout << "Entering C++ Script" << endl;

    std::random_device rd;
    rn gen(rd());

    GasContainer *box = new GasContainer(30, 1., 0.01, .01);
    box->evolve(100.);

    std::vector<std::vector<vd> > pos_animation = box->pos_animation;
    const int n1 = pos_animation.size();
    const int n2 = pos_animation[0].size();
    const int n3 = pos_animation[0][0].size();

    vd buffer;
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            for (int k = 0; k < n3; k++) {
                buffer.push_back(pos_animation[i][j][k]);
            }
        }
    }

    // hid_t file_id = H5Fcreate("pos_animation.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // hsize_t dims[3] = {n1, n2, n3};
    // hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

    // // Create a new dataset in the file
    // hid_t dataset_id = H5Dcreate2(file_id, "/pos_animation", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // // Write the data to the dataset
    // H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);    

    // // Close the dataset, dataspace, and file
    // H5Dclose(dataset_id);
    // H5Sclose(dataspace_id);
    // H5Fclose(file_id);

    hid_t file_id = H5Fcreate("pos_animation.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dims[3] = {static_cast<hsize_t>(n1), static_cast<hsize_t>(n2), static_cast<hsize_t>(n3)};
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);

    // Define a suitable datatype for your dataset (double)
    hid_t datatype_id = H5Tcopy(H5T_NATIVE_DOUBLE);

    // Create a new dataset in the file
    hid_t dataset_id = H5Dcreate2(file_id, "/pos_animation", datatype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Write the data to the dataset
    H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]);    

    // Close the dataset, dataspace, datatype, and file
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(datatype_id);
    H5Fclose(file_id);

    delete box;

    cout << "Exiting C++ Script" << endl;

    return 0;
}