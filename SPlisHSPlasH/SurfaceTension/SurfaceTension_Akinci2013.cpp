#include "SurfaceTension_Akinci2013.h"
#include <iostream>
#include "../Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;

namespace
{
	void addProjectedCapillaryAcceleration(
		Vector3r &ai,
		const Vector3r &boundaryNormal,
		const Real adhesionScale,
		const Real kernelWeight,
		const Real capillaryStrength,
		const Vector3r &capillaryDirection,
		const Real capillaryForwardStrength,
		const Vector3r &capillaryForwardDirection)
	{
		if ((kernelWeight <= static_cast<Real>(0.0)) || (adhesionScale <= static_cast<Real>(0.0)))
			return;

		const Real normalLength2 = boundaryNormal.squaredNorm();
		if (normalLength2 <= static_cast<Real>(1.0e-9))
			return;
		const Vector3r n = (static_cast<Real>(1.0) / sqrt(normalLength2)) * boundaryNormal;

		const auto addProjectedDirection = [&](const Vector3r &direction, const Real strength)
		{
			if (strength <= static_cast<Real>(0.0))
				return;

			const Real directionLength2 = direction.squaredNorm();
			if (directionLength2 <= static_cast<Real>(1.0e-9))
				return;

			Vector3r tangent = direction - direction.dot(n) * n;
			const Real tangentLength2 = tangent.squaredNorm();
			if (tangentLength2 <= static_cast<Real>(1.0e-9))
				return;

			tangent = (static_cast<Real>(1.0) / sqrt(tangentLength2)) * tangent;
			ai += adhesionScale * kernelWeight * strength * tangent;
		};

		addProjectedDirection(capillaryDirection, capillaryStrength);
		addProjectedDirection(capillaryForwardDirection, capillaryForwardStrength);
	}
}

SurfaceTension_Akinci2013::SurfaceTension_Akinci2013(FluidModel *model) :
	SurfaceTensionBase(model)
{
	m_normals.resize(model->numParticles(), Vector3r::Zero());

	model->addField({ "normal", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_normals[i][0]; }, false });
}

SurfaceTension_Akinci2013::~SurfaceTension_Akinci2013(void)
{
	m_model->removeFieldByName("normal");
	m_normals.clear();
}


void SurfaceTension_Akinci2013::computeNormals()
{
	Simulation *sim = Simulation::getCurrent();
	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel *model = m_model;

	// Compute normals
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r &ni = getNormal(i);
			ni.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = m_model->getDensity(neighborIndex);
				ni += m_model->getMass(neighborIndex) / density_j * sim->gradW(xi - xj);
			)
			ni = supportRadius*ni;
		}
	}

}

void SurfaceTension_Akinci2013::step()
{
	Simulation *sim = Simulation::getCurrent();
	const Real density0 = m_model->getDensity0();
	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	const Real kb = m_surfaceTensionBoundary; // Adhesion strength coefficient.
	const Real boundaryRepulsionDistance = static_cast<Real>(1.5) * sim->getParticleRadius();
	const Real boundaryRepulsionDistance2 = boundaryRepulsionDistance * boundaryRepulsionDistance;
	const Real boundaryRepulsionInvDistance = static_cast<Real>(1.0) / boundaryRepulsionDistance;
	const Real timeStepSize = TimeManager::getCurrent()->getTimeStepSize();
	const Real timeStepSize2 = timeStepSize * timeStepSize;
	const Real boundaryRepulsionScale = (timeStepSize2 > static_cast<Real>(1.0e-12)) ? static_cast<Real>(0.05) / timeStepSize2 : static_cast<Real>(0.0);
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	FluidModel *model = m_model;

	computeNormals();

	// Compute forces
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &ni = getNormal(i);
			const Real &rhoi = m_model->getDensity(i);
			Vector3r &ai = m_model->getAcceleration(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real &rhoj = m_model->getDensity(neighborIndex);
				const Real K_ij = static_cast<Real>(2.0)*density0 / (rhoi + rhoj);

				Vector3r accel;
				accel.setZero();

				// Cohesion force
				Vector3r xixj = (xi - xj);
				const Real length2 = xixj.squaredNorm();
				if (length2 > 1.0e-9)
				{
					xixj = (static_cast<Real>(1.0) / sqrt(length2)) * xixj;
					accel -= k * m_model->getMass(neighborIndex) * xixj * CohesionKernel::W(xi - xj);
				}

				// Curvature
				const Vector3r &nj = getNormal(neighborIndex);
				accel -= k * (ni - nj);

				ai += K_ij * accel;
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					if (bm_neighbor->isWall())
					{
						const Real boundaryAdhesionScale = bm_neighbor->getAdhesionScale();
						const Real boundaryAdhesion = kb * boundaryAdhesionScale;

						// adhesion force
						Vector3r xixj = (xi - xj);
						const Real length2 = xixj.squaredNorm();
						if (length2 > 1.0e-9)
						{
							const Real length = sqrt(length2);
							xixj = ((Real) 1.0 / length) * xixj;

							// Numerical guard against fluid particles entering particle boundaries.
							if (length2 < boundaryRepulsionDistance2)
							{
								const Real overlap = boundaryRepulsionDistance - length;
								const Real q = overlap * boundaryRepulsionInvDistance;
								const Real volumeRatio = bm_neighbor->getVolume(neighborIndex) / m_model->getVolume(i);
								ai += boundaryRepulsionScale * overlap * q * q * volumeRatio * xixj;
							}

							ai -= boundaryAdhesion * density0 * bm_neighbor->getVolume(neighborIndex) * xixj * AdhesionKernel::W(xi - xj);

							const Real capillaryWeight = (bm_neighbor->getVolume(neighborIndex) / m_model->getVolume(i)) * sim->W(xi - xj) / sim->W_zero();
							addProjectedCapillaryAcceleration(
								ai,
								xixj,
								boundaryAdhesionScale,
								capillaryWeight,
								bm_neighbor->getCapillaryStrength(),
								bm_neighbor->getCapillaryDirection(),
								bm_neighbor->getCapillaryForwardStrength(),
								bm_neighbor->getCapillaryForwardDirection());
						}
					}
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					if (bm_neighbor->isWall())
					{
						const Real boundaryAdhesionScale = bm_neighbor->getAdhesionScale();
						const Real boundaryAdhesion = kb * boundaryAdhesionScale;
						Vector3r xixj = xi - xj;
						const Real length2 = xixj.squaredNorm();
						if (length2 > 1.0e-9)
						{
							xixj = ((Real) 1.0 / sqrt(length2)) * xixj;
							ai -= boundaryAdhesion * density0 * xixj * rho * AdhesionKernel::W_zero() / sim->W_zero();

							const Real capillaryWeight = rho / density0;
							addProjectedCapillaryAcceleration(
								ai,
								xixj,
								boundaryAdhesionScale,
								capillaryWeight,
								bm_neighbor->getCapillaryStrength(),
								bm_neighbor->getCapillaryDirection(),
								bm_neighbor->getCapillaryForwardStrength(),
								bm_neighbor->getCapillaryForwardDirection());
						}
					}
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					if (bm_neighbor->isWall())
					{
						const Real boundaryAdhesionScale = bm_neighbor->getAdhesionScale();
						const Real boundaryAdhesion = kb * boundaryAdhesionScale;
						Vector3r xixj = (xi - xj);
						const Real length2 = xixj.squaredNorm();
						if (length2 > 1.0e-9)
						{
							xixj = ((Real) 1.0 / sqrt(length2)) * xixj;
							ai -= boundaryAdhesion * Vj * density0 * xixj * AdhesionKernel::W(xi - xj);

							const Real capillaryWeight = (Vj / m_model->getVolume(i)) * sim->W(xi - xj) / sim->W_zero();
							addProjectedCapillaryAcceleration(
								ai,
								xixj,
								boundaryAdhesionScale,
								capillaryWeight,
								bm_neighbor->getCapillaryStrength(),
								bm_neighbor->getCapillaryDirection(),
								bm_neighbor->getCapillaryForwardStrength(),
								bm_neighbor->getCapillaryForwardDirection());
						}
					}
				);
			}
		}
	}
}


void SurfaceTension_Akinci2013::reset()
{
}

void SurfaceTension_Akinci2013::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_normals[0]);
}

