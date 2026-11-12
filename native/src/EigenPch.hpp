#pragma once

// Eigen is deliberately confined to implementation targets.  Keeping this
// stable precompiled header free of DynSysLab headers means ordinary API and
// workflow changes do not invalidate the expensive template parse.
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
