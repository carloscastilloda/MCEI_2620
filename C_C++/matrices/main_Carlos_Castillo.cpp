// ==========================================================================
//  Taller de Algebra Lineal Numerica - MCEI
// ==========================================================================
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Reloj = std::chrono::high_resolution_clock;

// ==========================================================================
//  PUNTO 1 (C++ / Eigen3): Inversa de Moore-Penrose
//  Calcula A^+ via SVD y verifica las 4 propiedades de Penrose.
//  A es 4x3 y DEFICIENTE DE RANGO (rank = 2, col3 = 2*col1 - col2),
//  por lo que (A'A)^-1 A' no existe y hay que ir por SVD con truncacion.
// Pseudoinversa de Moore-Penrose via SVD con truncacion de valores singulares
//   A = U S V'   =>   A^+ = V S^+ U'
//   S^+ = diag(1/sigma_i) solo para sigma_i > tol; el resto se anula.
// --------------------------------------------------------------------------
MatrixXd pseudoInversa(const MatrixXd& A, double& tol_out, int& rango_out)
{
    Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const VectorXd& s = svd.singularValues();

    // Misma tolerancia que usa pinv() de Octave/MATLAB
    const double tol = std::max(A.rows(), A.cols()) *
                       s(0) * std::numeric_limits<double>::epsilon();
    tol_out = tol;

    VectorXd sinv = VectorXd::Zero(s.size());
    int r = 0;
    for (int i = 0; i < s.size(); ++i)
        if (s(i) > tol) { sinv(i) = 1.0 / s(i); ++r; }
    rango_out = r;

    return svd.matrixV() * sinv.asDiagonal() * svd.matrixU().transpose();
}

static const char* veredicto(double e) { return (e < 1e-10) ? "CUMPLE" : "NO CUMPLE"; }

int main()
{
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\n===== PUNTO 1 (C++/Eigen3): PSEUDOINVERSA DE MOORE-PENROSE =====\n";

    MatrixXd A(4, 3);
    A <<  1,  0,  2,
          2, -1,  5,
          0,  1, -1,
          1,  3, -1;

    std::cout << "\nA (" << A.rows() << "x" << A.cols() << ") =\n"
              << std::fixed << std::setprecision(4) << A << "\n";

    // SVD y rango
    Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const VectorXd& s = svd.singularValues();

    std::cout << std::scientific << std::setprecision(12);
    std::cout << "\nValores singulares:\n";
    for (int i = 0; i < s.size(); ++i)
        std::cout << "   sigma_" << i + 1 << " = " << s(i) << "\n";

    double tol; int r;
    MatrixXd Ap = pseudoInversa(A, tol, r);

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Tolerancia de rango  : " << tol << "\n";
    std::cout << "Rango numerico       : " << r << " (max " << std::min(A.rows(), A.cols()) << ")\n";
    std::cout << "cond_2(A) efectivo   : " << s(0) / s(r - 1) << "\n";
    if (r < std::min(A.rows(), A.cols())) {
        std::cout << ">> A es DEFICIENTE DE RANGO: A'A es singular,\n"
                  << "   la formula (A'A)^-1 A' NO es aplicable. Se usa SVD.\n";
    }

    std::cout << "\nA^+ (via SVD) =\n" << std::fixed << std::setprecision(8) << Ap << "\n";

    // comparacion con Eigen 3.3+
#if EIGEN_VERSION_AT_LEAST(3,3,0)
    MatrixXd Ap_cod = A.completeOrthogonalDecomposition().pseudoInverse();
    std::cout << std::scientific << std::setprecision(6)
              << "||A^+_SVD - A^+_COD(Eigen)||_F = "
              << (Ap - Ap_cod).norm() << "\n";
#endif

    // que pasa con las ecuaciones normales
    MatrixXd AtA = A.transpose() * A;
    std::cout << "\n--- Formula (A'A)^-1 A' ---\n";
    std::cout << "det(A'A) = " << AtA.determinant() << "   (~0 => singular)\n";
    Eigen::JacobiSVD<MatrixXd> svdN(AtA);
    std::cout << "cond_2(A'A) = "
              << svdN.singularValues()(0) /
                 svdN.singularValues()(svdN.singularValues().size() - 1) << "\n";
    MatrixXd Ap_ne = AtA.inverse() * A.transpose();
    std::cout << "||(A'A)^-1 A' - A^+||_F = " << (Ap_ne - Ap).norm()
              << "  <-- basura numerica\n";

    // 4 propiedades de Penrose
    const double e1 = (A * Ap * A - A).norm();
    const double e2 = (Ap * A * Ap - Ap).norm();
    const double e3 = ((A * Ap).transpose() - (A * Ap)).norm();
    const double e4 = ((Ap * A).transpose() - (Ap * A)).norm();

    std::cout << "\n--- VERIFICACION DE LAS 4 PROPIEDADES DE PENROSE ---\n";
    std::cout << "(1) ||A A+ A - A||_F     = " << e1 << "   " << veredicto(e1) << "\n";
    std::cout << "(2) ||A+ A A+ - A+||_F   = " << e2 << "   " << veredicto(e2) << "\n";
    std::cout << "(3) ||(A A+)' - A A+||_F = " << e3 << "   " << veredicto(e3) << "\n";
    std::cout << "(4) ||(A+ A)' - A+ A||_F = " << e4 << "   " << veredicto(e4) << "\n";

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\nA A^+ (proyector sobre Im(A)):\n" << A * Ap
              << "\ntraza = " << (A * Ap).trace() << "  (= rango)\n";
    std::cout << "\nA^+ A (proyector sobre Im(A')):\n" << Ap * A
              << "\ntraza = " << (Ap * A).trace() << "  (= rango)\n";
    std::cout << "\nNota: A^+A != I porque rank(A) = " << r << " < n = " << A.cols() << ".\n";
    return 0;
}
// ==========================================================================
//  PUNTO 2 (C++ / Eigen3): sistema lineal 10x10 por Gauss, LU y QR
//
//  El enunciado da A pero no b. Se toma b = A*ones(10,1) para conocer la
//  solucion exacta (todo unos) y poder medir el ERROR, no solo el residual.
// ==========================================================================
// --------------------------------------------------------------------------
//  Eliminacion de Gauss con pivoteo parcial, implementada explicitamente
// --------------------------------------------------------------------------
VectorXd gaussPP(MatrixXd A, VectorXd b, int& npiv)
{
    const int n = static_cast<int>(A.rows());
    npiv = 0;

    for (int k = 0; k < n - 1; ++k) {
        // pivoteo parcial: mayor |elemento| de la columna k desde la fila k
        int p = k;
        double maxv = std::abs(A(k, k));
        for (int i = k + 1; i < n; ++i)
            if (std::abs(A(i, k)) > maxv) { maxv = std::abs(A(i, k)); p = i; }

        if (p != k) {
            A.row(k).swap(A.row(p));
            std::swap(b(k), b(p));
            ++npiv;
        }
        if (std::abs(A(k, k)) < std::numeric_limits<double>::epsilon())
            throw std::runtime_error("Matriz singular en la columna " + std::to_string(k));

        for (int i = k + 1; i < n; ++i) {
            const double f = A(i, k) / A(k, k);
            A.row(i).tail(n - k) -= f * A.row(k).tail(n - k);
            b(i) -= f * b(k);
            A(i, k) = 0.0;
        }
    }

    // sustitucion regresiva
    VectorXd x(n);
    for (int i = n - 1; i >= 0; --i) {
        double suma = b(i);
        for (int j = i + 1; j < n; ++j) suma -= A(i, j) * x(j);
        x(i) = suma / A(i, i);
    }
    return x;
}

// cond_2 via SVD 
double cond2(const MatrixXd& M)
{
    Eigen::JacobiSVD<MatrixXd> svd(M);
    const auto& s = svd.singularValues();
    return s(0) / s(s.size() - 1);
}

int main()
{
    std::cout << "\n===== PUNTO 2 (C++/Eigen3): SISTEMA LINEAL 10x10 =====\n";

    MatrixXd A(10, 10);
    A << 2, 1, 0, 3, 2, 1, 0, 2, 1, 4,
         1, 3, 2, 0, 1, 4, 2, 1, 0, 2,
         0, 2, 4, 1, 3, 0, 1, 2, 4, 1,
         3, 0, 1, 5, 2, 1, 3, 0, 2, 1,
         2, 1, 3, 2, 6, 2, 1, 4, 0, 3,
         1, 4, 0, 1, 2, 5, 2, 1, 3, 0,
         0, 2, 1, 3, 1, 2, 4, 0, 2, 1,
         2, 1, 2, 0, 4, 1, 0, 5, 3, 2,
         1, 0, 4, 2, 0, 3, 2, 3, 6, 1,
         4, 2, 1, 1, 3, 0, 1, 2, 1, 5;

    const int n = 10;
    VectorXd x_exact = VectorXd::Ones(n);
    VectorXd b       = A * x_exact;

    std::cout << std::fixed << std::setprecision(0) << "b' = " << b.transpose() << "\n";

    // caracterizacion de A
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n--- Caracterizacion de A ---\n";
    std::cout << "det(A)      = " << A.determinant() << "\n";
    std::cout << "cond_2(A)   = " << cond2(A) << "\n";
    std::cout << "Simetrica?  : " << (A.isApprox(A.transpose()) ? "si" : "no") << "\n";
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(A);
    std::cout << "Autovalores : " << es.eigenvalues().transpose() << "\n";
    if (es.eigenvalues().minCoeff() > 0)
        std::cout << "A es definida positiva -> Cholesky aplicable.\n";
    else
        std::cout << "A es simetrica pero INDEFINIDA: Cholesky NO aplicable.\n";

    // los tres metodos
    int npiv = 0;
    auto t0 = Reloj::now();
    VectorXd x_gauss = gaussPP(A, b, npiv);
    double t_gauss = std::chrono::duration<double>(Reloj::now() - t0).count();

    t0 = Reloj::now();
    Eigen::PartialPivLU<MatrixXd> lu(A);       // PA = LU
    VectorXd x_lu = lu.solve(b);
    double t_lu = std::chrono::duration<double>(Reloj::now() - t0).count();

    t0 = Reloj::now();
    Eigen::HouseholderQR<MatrixXd> qr(A);      // A = QR
    VectorXd x_qr = qr.solve(b);
    double t_qr = std::chrono::duration<double>(Reloj::now() - t0).count();

    t0 = Reloj::now();
    Eigen::FullPivLU<MatrixXd> flu(A);         // referencia robusta
    VectorXd x_ref = flu.solve(b);
    double t_ref = std::chrono::duration<double>(Reloj::now() - t0).count();

    // resultados
    std::cout << "\n--- Soluciones (deben ser todas 1) ---\n";
    std::cout << "  i |    Gauss    |     LU      |     QR      |  FullPivLU\n";
    for (int i = 0; i < n; ++i)
        std::cout << std::setw(3) << i + 1 << " | "
                  << std::setw(11) << std::setprecision(8) << x_gauss(i) << " | "
                  << std::setw(11) << x_lu(i) << " | "
                  << std::setw(11) << x_qr(i) << " | "
                  << std::setw(11) << x_ref(i) << "\n";

    std::cout << "\n--- Comparacion de estabilidad y eficiencia ---\n";
    std::cout << std::left << std::setw(14) << "Metodo"
              << std::setw(16) << "||Ax-b||_2" << std::setw(16) << "err.rel."
              << std::setw(16) << "tiempo [s]" << "costo\n";

    const char* nom[4]   = {"Gauss(mano)", "PartialPivLU", "HouseholderQR", "FullPivLU"};
    const char* costo[4] = {"2n^3/3", "2n^3/3", "4n^3/3", "2n^3/3 + pivoteo"};
    VectorXd X[4] = {x_gauss, x_lu, x_qr, x_ref};
    double  tt[4] = {t_gauss, t_lu, t_qr, t_ref};
    std::cout << std::scientific << std::setprecision(4);
    for (int k = 0; k < 4; ++k)
        std::cout << std::left << std::setw(14) << nom[k]
                  << std::setw(16) << (A * X[k] - b).norm()
                  << std::setw(16) << (X[k] - x_exact).norm() / x_exact.norm()
                  << std::setw(16) << tt[k] << costo[k] << "\n";

    std::cout << "\nCota teorica: err_rel <= cond(A)*eps = "
              << cond2(A) * std::numeric_limits<double>::epsilon() << "\n";

    // factores obtenidos
    std::cout << "\n--- Factores ---\n";
    std::cout << "Permutaciones de fila en Gauss: " << npiv << "\n";
    MatrixXd LU = lu.matrixLU();
    MatrixXd L  = MatrixXd::Identity(n, n);
    L.triangularView<Eigen::StrictlyLower>() = LU.triangularView<Eigen::StrictlyLower>();
    MatrixXd U  = LU.triangularView<Eigen::Upper>();
    MatrixXd P  = lu.permutationP();
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "diag(U) = " << U.diagonal().transpose() << "\n";
    MatrixXd Rm = qr.matrixQR().triangularView<Eigen::Upper>();
    MatrixXd Qm = qr.householderQ();
    std::cout << "diag(R) = " << Rm.diagonal().transpose() << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "||P*A - L*U||_F = " << (P * A - L * U).norm() << "\n";
    std::cout << "||A - Q*R||_F   = " << (A - Qm * Rm).norm() << "\n";
    std::cout << "||Q'Q - I||_F   = "
              << (Qm.transpose() * Qm - MatrixXd::Identity(n, n)).norm() << "\n";

    // eficiencia promediada
    const int REP = 20000;
    std::cout << "\n--- Tiempo medio sobre " << REP << " repeticiones ---\n";
    double acc = 0;
    t0 = Reloj::now();
    for (int k = 0; k < REP; ++k) { int np; acc += gaussPP(A, b, np)(0); }
    double tg = std::chrono::duration<double>(Reloj::now() - t0).count() / REP;

    t0 = Reloj::now();
    for (int k = 0; k < REP; ++k) acc += Eigen::PartialPivLU<MatrixXd>(A).solve(b)(0);
    double tl = std::chrono::duration<double>(Reloj::now() - t0).count() / REP;

    t0 = Reloj::now();
    for (int k = 0; k < REP; ++k) acc += Eigen::HouseholderQR<MatrixXd>(A).solve(b)(0);
    double tq = std::chrono::duration<double>(Reloj::now() - t0).count() / REP;

    std::cout << "Gauss (a mano)  : " << tg << " s\n";
    std::cout << "PartialPivLU    : " << tl << " s\n";
    std::cout << "HouseholderQR   : " << tq << " s\n";
    std::cout << "(acumulador anti-optimizacion: " << acc << ")\n";

    std::cout << "\nVentaja de LU: con varios b se factoriza UNA vez (O(n^3))\n"
                 "y cada nuevo b cuesta solo O(n^2):\n";
    VectorXd b2 = VectorXd::LinSpaced(n, 1, n);
    std::cout << "  x(b2) con la MISMA factorizacion = " << std::fixed
              << std::setprecision(4) << lu.solve(b2).transpose() << "\n";
    return 0;
}

// ==========================================================================
//  PUNTO 3 (C++ / Eigen3): condicionamiento numerico y estabilidad
//   - cond() de una matriz mal condicionada (Hilbert)
//   - solucion con la matriz original y con una perturbacion pequena
//   - analisis de la amplificacion del error
// ==========================================================================
MatrixXd hilbert(int n)
{
    MatrixXd H(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            H(i, j) = 1.0 / (i + j + 1);   // H(i,j) = 1/(i+j-1) en base 1
    return H;
}

double cond2(const MatrixXd& M)
{
    Eigen::JacobiSVD<MatrixXd> svd(M);
    const auto& s = svd.singularValues();
    return s(0) / s(s.size() - 1);
}

int main()
{
    const double EPS = std::numeric_limits<double>::epsilon();
    std::cout << "\n===== PUNTO 3 (C++/Eigen3): CONDICIONAMIENTO Y ESTABILIDAD =====\n";

    //  3.1 numero de condicion
    std::cout << "\n--- 3.1 cond_2 de la matriz de Hilbert ---\n";
    std::cout << "  n |      cond_2(H)      | digitos perdidos | utiles (de 16)\n";
    for (int n = 2; n <= 12; ++n) {
        double c = cond2(hilbert(n));
        std::cout << std::setw(3) << n << " | " << std::scientific << std::setprecision(4)
                  << std::setw(19) << c << " | " << std::fixed << std::setprecision(1)
                  << std::setw(16) << std::log10(c) << " | "
                  << std::setw(6) << 16 - std::log10(c) << "\n";
    }

    const int n = 10;
    MatrixXd H = hilbert(n);
    Eigen::JacobiSVD<MatrixXd> svdH(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const VectorXd s = svdH.singularValues();

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "\nMatriz de estudio: H = hilbert(" << n << ")\n";
    std::cout << "cond_2(H) = " << cond2(H) << "\n";
    std::cout << "Valores singulares:\n";
    for (int i = 0; i < s.size(); ++i) std::cout << "   " << s(i) << "\n";

    //  3.2 sistema original vs perturbado
    VectorXd x_exact = VectorXd::Ones(n);
    VectorXd b       = H * x_exact;

    const double delta = 1e-10;
    VectorXd db     = VectorXd::Constant(n, delta);
    VectorXd b_pert = b + db;

    // se resuelve con LU con pivoteo completo (lo mas robusto disponible)
    Eigen::FullPivLU<MatrixXd> lu(H);
    VectorXd x0 = lu.solve(b);
    VectorXd x1 = lu.solve(b_pert);

    const double rel_b = db.norm() / b.norm();
    const double rel_x = (x1 - x0).norm() / x0.norm();

    std::cout << "\n--- 3.2 Perturbacion del lado derecho b ---\n";
    std::cout << "Perturbacion RELATIVA en b   : " << rel_b << "\n";
    std::cout << "Cambio RELATIVO en x         : " << rel_x << "\n";
    std::cout << "Factor de amplificacion real : " << rel_x / rel_b << "\n";
    std::cout << "Cota teorica cond_2(H)       : " << cond2(H) << "\n";

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "\n  i |   x original   |  x perturbado  |   diferencia\n";
    for (int i = 0; i < n; ++i)
        std::cout << std::setw(3) << i + 1 << " | " << std::setw(14) << x0(i)
                  << " | " << std::setw(14) << x1(i) << " | "
                  << std::scientific << std::setprecision(6) << std::setw(14)
                  << x1(i) - x0(i) << std::fixed << std::setprecision(8) << "\n";

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\nError relativo vs x_exacta (sin perturbar) : "
              << (x0 - x_exact).norm() / x_exact.norm() << "\n";
    std::cout << "Error relativo vs x_exacta (perturbado)    : "
              << (x1 - x_exact).norm() / x_exact.norm() << "\n";
    std::cout << "NOTA: incluso SIN perturbar hay error: almacenar H en doble\n"
                 "      precision ya introduce una perturbacion de ~1e-16.\n";

    //  contraste con una matriz BIEN condicionada
    MatrixXd Ab(10, 10);
    Ab << 2,1,0,3,2,1,0,2,1,4,  1,3,2,0,1,4,2,1,0,2,  0,2,4,1,3,0,1,2,4,1,
          3,0,1,5,2,1,3,0,2,1,  2,1,3,2,6,2,1,4,0,3,  1,4,0,1,2,5,2,1,3,0,
          0,2,1,3,1,2,4,0,2,1,  2,1,2,0,4,1,0,5,3,2,  1,0,4,2,0,3,2,3,6,1,
          4,2,1,1,3,0,1,2,1,5;
    VectorXd b_ok  = Ab * x_exact;
    VectorXd x0_ok = Ab.partialPivLu().solve(b_ok);
    VectorXd x1_ok = Ab.partialPivLu().solve(b_ok + db);
    const double rb = db.norm() / b_ok.norm();
    const double rx = (x1_ok - x0_ok).norm() / x0_ok.norm();
    std::cout << "\n--- Contraste con la matriz BIEN condicionada (punto 2) ---\n";
    std::cout << "cond_2 = " << cond2(Ab) << " | pert.rel.b = " << rb
              << " -> cambio rel.x = " << rx << " (amplif. " << rx / rb << ")\n";
    std::cout << "La MISMA perturbacion es inofensiva: el condicionamiento es\n"
                 "una propiedad del PROBLEMA, no del algoritmo.\n";

    //  3.3 estabilidad de los tres metodos
    std::cout << "\n--- 3.3 Estabilidad de los metodos ante el mal condicionamiento ---\n";
    VectorXd x_lu   = H.partialPivLu().solve(b);
    VectorXd x_lu2  = H.partialPivLu().solve(b_pert);
    VectorXd x_qr   = H.householderQr().solve(b);
    VectorXd x_qr2  = H.householderQr().solve(b_pert);
    VectorXd x_svd  = svdH.solve(b);
    VectorXd x_svd2 = svdH.solve(b_pert);

    const char* nom[3] = {"LU (Gauss)", "QR", "SVD"};
    VectorXd Xa[3] = {x_lu, x_qr, x_svd};
    VectorXd Xb[3] = {x_lu2, x_qr2, x_svd2};
    std::cout << std::left << std::setw(14) << "Metodo" << std::setw(16) << "err.rel.x"
              << std::setw(16) << "||Hx-b||" << "sensib.pert.\n";
    for (int j = 0; j < 3; ++j)
        std::cout << std::left << std::setw(14) << nom[j]
                  << std::setw(16) << (Xa[j] - x_exact).norm() / x_exact.norm()
                  << std::setw(16) << (H * Xa[j] - b).norm()
                  << (Xb[j] - Xa[j]).norm() / Xa[j].norm() << "\n";

    // SVD truncada = regularizacion
    const double tolr = 1e-12 * s(0);
    int k = 0; for (int i = 0; i < s.size(); ++i) if (s(i) > tolr) ++k;
    VectorXd Uty = svdH.matrixU().leftCols(k).transpose() * b;
    VectorXd x_tr = svdH.matrixV().leftCols(k) *
                    (Uty.array() / s.head(k).array()).matrix();
    std::cout << std::left << std::setw(14) << "SVD truncada"
              << std::setw(16) << (x_tr - x_exact).norm() / x_exact.norm()
              << std::setw(16) << (H * x_tr - b).norm()
              << "(rango efectivo " << k << "/" << n << ")\n";

    std::cout << "\nLECTURA:\n"
                 " * El RESIDUAL es minusculo en los tres metodos: son estables\n"
                 "   hacia atras (backward stable).\n"
                 " * El ERROR ||x - x_exacta|| es grande en los tres: el PROBLEMA\n"
                 "   esta mal condicionado y ningun algoritmo lo puede arreglar.\n"
                 " * Residual pequeno NO implica solucion correcta si cond(A) es alto.\n";
    std::cout << "\neps de la maquina = " << EPS << "\n";
    return 0;
}

// ==========================================================================
//  Taller de Algebra Lineal Numerica - MCEI
//  PUNTO 4 (C++ / Eigen3): inversa por metodo directo, QR y SVD
// ==========================================================================
MatrixXd hilbert(int n)
{
    MatrixXd H(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j) H(i, j) = 1.0 / (i + j + 1);
    return H;
}

double cond2(const MatrixXd& M)
{
    Eigen::JacobiSVD<MatrixXd> svd(M);
    const auto& s = svd.singularValues();
    return s(0) / s(s.size() - 1);
}

int main()
{
    const double EPS = std::numeric_limits<double>::epsilon();
    std::cout << "\n===== PUNTO 4 (C++/Eigen3): INVERSA DIRECTA, QR Y SVD =====\n";

    MatrixXd M1(10, 10);
    M1 << 2,1,0,3,2,1,0,2,1,4,  1,3,2,0,1,4,2,1,0,2,  0,2,4,1,3,0,1,2,4,1,
          3,0,1,5,2,1,3,0,2,1,  2,1,3,2,6,2,1,4,0,3,  1,4,0,1,2,5,2,1,3,0,
          0,2,1,3,1,2,4,0,2,1,  2,1,2,0,4,1,0,5,3,2,  1,0,4,2,0,3,2,3,6,1,
          4,2,1,1,3,0,1,2,1,5;
    MatrixXd M2 = hilbert(10);

    MatrixXd casos[2]      = {M1, M2};
    const char* nombres[2] = {"A del punto 2 (bien condicionada)",
                              "Hilbert 10x10 (mal condicionada)"};

    for (int c = 0; c < 2; ++c) {
        const MatrixXd& A = casos[c];
        const int n = static_cast<int>(A.rows());
        const MatrixXd I = MatrixXd::Identity(n, n);

        std::cout << "CASO " << c + 1 << ": " << nombres[c] << "\n";
        std::cout << std::scientific << std::setprecision(6)
                  << "cond_2(A) = " << cond2(A) << "\n";

        //  (a) inversa directa
        auto t0 = Reloj::now();
        MatrixXd Ainv_dir = A.inverse();
        double t_dir = std::chrono::duration<double>(Reloj::now() - t0).count();

        //  (b) inversa via QR
        // A = QR  =>  A^-1 = R^-1 Q'. Se resuelve R*Ainv = Q' por sustitucion
        // regresiva; nunca se forma R^-1 de manera explicita.
        t0 = Reloj::now();
        Eigen::HouseholderQR<MatrixXd> qr(A);
        MatrixXd Qt = qr.householderQ().transpose();
        MatrixXd R  = qr.matrixQR().triangularView<Eigen::Upper>();
        MatrixXd Ainv_qr = R.triangularView<Eigen::Upper>().solve(Qt);
        double t_qr = std::chrono::duration<double>(Reloj::now() - t0).count();

        //  (c) inversa via SVD
        t0 = Reloj::now();
        Eigen::JacobiSVD<MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const VectorXd s = svd.singularValues();
        MatrixXd Ainv_svd = svd.matrixV() *
                            s.cwiseInverse().asDiagonal() *
                            svd.matrixU().transpose();
        double t_svd = std::chrono::duration<double>(Reloj::now() - t0).count();

        // SVD truncada (pseudoinversa regularizada) con dos tolerancias:
        //   tol_std = tolerancia estandar de pinv()  -> normalmente NO trunca
        //   tol_agr = tolerancia agresiva            -> regulariza de verdad
        auto truncar = [&](double tol, int& k) {
            k = 0;
            for (int i = 0; i < s.size(); ++i) if (s(i) > tol) ++k;
            return MatrixXd(svd.matrixV().leftCols(k) *
                            s.head(k).cwiseInverse().asDiagonal() *
                            svd.matrixU().leftCols(k).transpose());
        };
        int k_std = 0, k_agr = 0;
        MatrixXd Ainv_tr  = truncar(n * s(0) * EPS, k_std);
        MatrixXd Ainv_tra = truncar(1e-10 * s(0),   k_agr);

        //  comparacion
        std::cout << "\n" << std::left << std::setw(18) << "Metodo"
                  << std::setw(16) << "||A*Ai-I||_F" << std::setw(16) << "||Ai*A-I||_F"
                  << std::setw(18) << "dif.rel vs inv()" << "tiempo[s]\n";
        const char* met[3] = {"Directo inverse", "QR", "SVD"};
        MatrixXd Ai[3] = {Ainv_dir, Ainv_qr, Ainv_svd};
        double tt[3]   = {t_dir, t_qr, t_svd};
        const double nrm = Ainv_dir.norm();
        std::cout << std::scientific << std::setprecision(4);
        for (int j = 0; j < 3; ++j)
            std::cout << std::left << std::setw(18) << met[j]
                      << std::setw(16) << (A * Ai[j] - I).norm()
                      << std::setw(16) << (Ai[j] * A - I).norm()
                      << std::setw(18) << (Ai[j] - Ainv_dir).norm() / nrm
                      << tt[j] << "\n";
        std::cout << std::left << std::setw(18) << "SVD trunc (std)"
                  << std::setw(16) << (A * Ainv_tr - I).norm()
                  << std::setw(16) << (Ainv_tr * A - I).norm()
                  << "(rango efectivo " << k_std << "/" << n << ")\n";
        std::cout << std::left << std::setw(18) << "SVD trunc (1e-10)"
                  << std::setw(16) << (A * Ainv_tra - I).norm()
                  << std::setw(16) << (Ainv_tra * A - I).norm()
                  << "(rango efectivo " << k_agr << "/" << n << ")\n";

        std::cout << "\nCota teorica del error relativo: cond_2(A)*eps = "
                  << cond2(A) * EPS << "\n";

        //  prueba de uso real: x = Ainv * b con b = A*ones 
        // ||A*Ai - I|| castiga a la version truncada (A*Ai es un proyector de
        // rango k, asi que su distancia a I es sqrt(n-k) por construccion).
        // La metrica que importa es el ERROR en la solucion del sistema.
        const VectorXd x_ex = VectorXd::Ones(n);
        const VectorXd bvec = A * x_ex;
        const char* metx[5] = {"Directo inverse", "QR", "SVD",
                               "SVD trunc (std)", "SVD trunc (1e-10)"};
        MatrixXd Aix[5] = {Ainv_dir, Ainv_qr, Ainv_svd, Ainv_tr, Ainv_tra};
        std::cout << "\nError al RESOLVER A x = b con cada inversa (b = A*ones):\n";
        for (int j = 0; j < 5; ++j)
            std::cout << "   " << std::left << std::setw(18) << metx[j]
                      << "err.rel(x) = " << (Aix[j] * bvec - x_ex).norm() / x_ex.norm()
                      << "\n";

        if (c == 1) {
            std::cout << "\nValores singulares de Hilbert(10):\n";
            for (int i = 0; i < s.size(); ++i) std::cout << "   " << s(i) << "\n";
            std::cout << "Los ultimos sigma estan al nivel del ruido de redondeo;\n"
                         "al invertirlos (1/sigma) ese ruido se amplifica ~1e13.\n";
        }
    }

    std::cout << "\n--- CONCLUSION DEL PUNTO 4 ---\n"
                 "1) Con matriz BIEN condicionada los tres metodos coinciden a\n"
                 "   nivel de ~1e-15; gana la inversa directa por costo (2n^3\n"
                 "   frente a ~4n^3/3+n^3 de QR y ~12-22n^3 de la SVD).\n"
                 "2) Con matriz MAL condicionada NINGUNO de los tres se salva:\n"
                 "   los tres errores quedan del orden de cond(A)*eps, que es la\n"
                 "   cota teorica. Esa perdida de precision es del PROBLEMA, no\n"
                 "   del algoritmo, y por eso cambiar de metodo no la elimina.\n"
                 "3) La ventaja real de la SVD no es un error mas pequeno sino que\n"
                 "   es el unico metodo que DIAGNOSTICA (muestra los sigma\n"
                 "   pequenos) y permite REGULARIZAR truncandolos: la ultima fila\n"
                 "   de la tabla del caso 2 lo hace explicito.\n"
                 "4) REGLA DE ORO: para resolver A x = b NUNCA se calcula A^-1.\n"
                 "   Un solve() es mas rapido y mas preciso.\n";
    return 0;
}

// ==========================================================================
//  Taller de Algebra Lineal Numerica - MCEI
//  PROBLEMA PRACTICO (C++ / Eigen3): convertidor DC-DC
//
//  P = b0 + b1*V + b2*I + b3*T + e      ->      y = X*beta
//
//  Minimos cuadrados por 4 vias: ecuaciones normales, QR, SVD y
//  BDCSVD/COD; residual, MSE/RMSE, R^2, influencia relativa y cond(X).
//
// ==========================================================================
// --------------------------------------------------------------------------
//  Lector de CSV con encabezado (devuelve la matriz de datos n x 4)
// --------------------------------------------------------------------------
MatrixXd leerCSV(const std::string& ruta)
{
    std::ifstream f(ruta);
    if (!f) throw std::runtime_error("No se pudo abrir el archivo: " + ruta);

    std::string linea;
    std::getline(f, linea);                       // descarta el encabezado

    std::vector<std::vector<double>> filas;
    while (std::getline(f, linea)) {
        if (linea.empty()) continue;
        std::stringstream ss(linea);
        std::string campo;
        std::vector<double> fila;
        while (std::getline(ss, campo, ',')) fila.push_back(std::stod(campo));
        if (!fila.empty()) filas.push_back(fila);
    }
    if (filas.empty()) throw std::runtime_error("CSV vacio");

    MatrixXd D(filas.size(), filas[0].size());
    for (size_t i = 0; i < filas.size(); ++i)
        for (size_t j = 0; j < filas[i].size(); ++j) D(i, j) = filas[i][j];
    return D;
}

double cond2(const MatrixXd& M)
{
    Eigen::JacobiSVD<MatrixXd> svd(M);
    const auto& s = svd.singularValues();
    return s(0) / s(s.size() - 1);
}

double desvEst(const VectorXd& v)                       // desviacion muestral
{
    const double mu = v.mean();
    return std::sqrt((v.array() - mu).square().sum() / (v.size() - 1));
}

double correl(const VectorXd& a, const VectorXd& b)
{
    VectorXd x = a.array() - a.mean();
    VectorXd y = b.array() - b.mean();
    return x.dot(y) / (x.norm() * y.norm());
}

int main(int argc, char** argv)
{
    const std::string ruta = (argc > 1) ? argv[1]
                                        : "/home/carlos/MCEI_2620/C_C++/matrices/datos_convertidor_realista.csv";
    std::cout << "\n===== PROBLEMA PRACTICO (C++/Eigen3): CONVERTIDOR DC-DC =====\n";

    MatrixXd D = leerCSV(ruta);
    const int n = static_cast<int>(D.rows());
    const int p = 4;

    VectorXd V  = D.col(0);      // voltaje [V]
    VectorXd Ic = D.col(1);      // corriente [A]
    VectorXd T  = D.col(2);      // temperatura [C]
    VectorXd y  = D.col(3);      // potencia [W]

    std::cout << "Archivo : " << ruta << "\n";
    std::cout << "Muestras: n = " << n << "\n";

    const char* etq[4] = {"Voltaje [V]", "Corriente [A]", "Temp [C]", "Potencia [W]"};
    std::cout << "\n" << std::left << std::setw(16) << "Variable" << std::right
              << std::setw(11) << "media" << std::setw(11) << "desv.est"
              << std::setw(11) << "min" << std::setw(11) << "max" << "\n";
    std::cout << std::fixed << std::setprecision(4);
    for (int k = 0; k < 4; ++k)
        std::cout << std::left << std::setw(16) << etq[k] << std::right
                  << std::setw(11) << D.col(k).mean()
                  << std::setw(11) << desvEst(D.col(k))
                  << std::setw(11) << D.col(k).minCoeff()
                  << std::setw(11) << D.col(k).maxCoeff() << "\n";

    //  matriz de diseno
    MatrixXd X(n, p);
    X.col(0) = VectorXd::Ones(n);
    X.col(1) = V;
    X.col(2) = Ic;
    X.col(3) = T;

    //  (a) ecuaciones normales
    MatrixXd XtX = X.transpose() * X;
    VectorXd Xty = X.transpose() * y;
    auto t0 = Reloj::now();
    VectorXd b_ne = XtX.ldlt().solve(Xty);            // XtX es simetrica def.pos.
    double t_ne = std::chrono::duration<double>(Reloj::now() - t0).count();
    VectorXd b_ne_inv = XtX.inverse() * Xty;          // version con inversa explicita

    //  (b) QR (Householder)
    t0 = Reloj::now();
    VectorXd b_qr = X.householderQr().solve(y);
    double t_qr = std::chrono::duration<double>(Reloj::now() - t0).count();

    //  (c) SVD
    t0 = Reloj::now();
    Eigen::JacobiSVD<MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);
    VectorXd b_svd = svd.solve(y);
    double t_svd = std::chrono::duration<double>(Reloj::now() - t0).count();

    //  (d) QR con pivoteo de columnas (robusta)
    t0 = Reloj::now();
    VectorXd b_cpqr = X.colPivHouseholderQr().solve(y);
    double t_cp = std::chrono::duration<double>(Reloj::now() - t0).count();

    //  ENTREGABLE 1: coeficientes
    std::cout << "\n===== 1. COEFICIENTES ESTIMADOS =====\n";
    const char* nomb[4] = {"b0 (interc)", "b1 (V)", "b2 (I)", "b3 (T)"};
    std::cout << std::left << std::setw(13) << "coef" << std::right
              << std::setw(15) << "Ec.normales" << std::setw(15) << "QR"
              << std::setw(15) << "SVD" << std::setw(15) << "ColPivQR" << "\n";
    std::cout << std::fixed << std::setprecision(8);
    for (int k = 0; k < p; ++k)
        std::cout << std::left << std::setw(13) << nomb[k] << std::right
                  << std::setw(15) << b_ne(k) << std::setw(15) << b_qr(k)
                  << std::setw(15) << b_svd(k) << std::setw(15) << b_cpqr(k) << "\n";

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "\nDiferencias respecto a QR (norma 2):\n";
    std::cout << "  Ec. normales (LDLT): " << (b_ne - b_qr).norm() << "\n";
    std::cout << "  Ec. normales inv() : " << (b_ne_inv - b_qr).norm() << "\n";
    std::cout << "  SVD                : " << (b_svd - b_qr).norm() << "\n";
    std::cout << "  ColPivQR           : " << (b_cpqr - b_qr).norm() << "\n";
    std::cout << "Tiempos [s]: normales " << t_ne << " | QR " << t_qr
              << " | SVD " << t_svd << " | ColPivQR " << t_cp << "\n";

    const VectorXd beta = b_qr;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\n>>> MODELO AJUSTADO:\n    P = " << beta(0) << std::showpos
              << " " << beta(1) << "*V " << beta(2) << "*I "
              << beta(3) << "*T" << std::noshowpos << "\n";

    //  ENTREGABLE 3: residual y error cuadratico
    VectorXd r   = y - X * beta;
    const double SSE  = r.squaredNorm();
    const double MSE  = SSE / n;
    const double RMSE = std::sqrt(MSE);
    const double MAE  = r.cwiseAbs().mean();
    const double SST  = (y.array() - y.mean()).square().sum();
    const double R2   = 1.0 - SSE / SST;
    const double R2a  = 1.0 - (1.0 - R2) * (n - 1.0) / (n - p);

    std::cout << "\n===== 3. RESIDUAL Y ERROR CUADRATICO MEDIO =====\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "||r||_2 = " << r.norm() << " W\n";
    std::cout << "SSE     = " << SSE  << " W^2\n";
    std::cout << "MSE     = " << MSE  << " W^2\n";
    std::cout << "RMSE    = " << RMSE << " W\n";
    std::cout << "MAE     = " << MAE  << " W\n";
    std::cout << std::setprecision(8);
    std::cout << "R^2          = " << R2  << "\n";
    std::cout << "R^2 ajustado = " << R2a << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "media(r) = " << r.mean() << "   (debe ser ~0)\n";
    std::cout << "||X'r||_inf = " << (X.transpose() * r).cwiseAbs().maxCoeff()
              << "   (condicion de ortogonalidad de MC)\n";
    std::cout << std::fixed << std::setprecision(4)
              << "max|r| = " << r.cwiseAbs().maxCoeff() << " W    RMSE/media(P) = "
              << 100.0 * RMSE / y.mean() << " %\n";

    //  ENTREGABLE 2: influencia relativa de cada variable
    const double sdV = desvEst(V), sdI = desvEst(Ic), sdT = desvEst(T), sdY = desvEst(y);
    const double bstd[3] = {beta(1) * sdV / sdY, beta(2) * sdI / sdY, beta(3) * sdT / sdY};
    const double suma = std::abs(bstd[0]) + std::abs(bstd[1]) + std::abs(bstd[2]);

    const double sigma2 = SSE / (n - p);
    MatrixXd C  = sigma2 * XtX.inverse();             // covarianza de beta
    VectorXd se = C.diagonal().cwiseSqrt();

    std::cout << "\n===== 2. INFLUENCIA RELATIVA DE CADA VARIABLE =====\n";
    std::cout << std::left << std::setw(12) << "var" << std::right
              << std::setw(13) << "beta" << std::setw(13) << "err.est"
              << std::setw(11) << "t" << std::setw(13) << "beta_std"
              << std::setw(13) << "influencia%" << "\n";
    std::cout << std::fixed << std::setprecision(6);
    for (int k = 0; k < 3; ++k)
        std::cout << std::left << std::setw(12) << nomb[k + 1] << std::right
                  << std::setw(13) << beta(k + 1) << std::setw(13) << se(k + 1)
                  << std::setw(11) << std::setprecision(4) << beta(k + 1) / se(k + 1)
                  << std::setw(13) << std::setprecision(6) << bstd[k]
                  << std::setw(13) << std::setprecision(2) << 100 * std::abs(bstd[k]) / suma
                  << std::setprecision(6) << "\n";
    std::cout << std::left << std::setw(12) << nomb[0] << std::right
              << std::setw(13) << beta(0) << std::setw(13) << se(0)
              << std::setw(11) << std::setprecision(4) << beta(0) / se(0) << "\n";

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nEfecto sobre P al recorrer TODO el rango de cada variable:\n";
    std::cout << "  V: " << std::showpos << beta(1) * (V.maxCoeff() - V.minCoeff())
              << std::noshowpos << " W\n";
    std::cout << "  I: " << std::showpos << beta(2) * (Ic.maxCoeff() - Ic.minCoeff())
              << std::noshowpos << " W\n";
    std::cout << "  T: " << std::showpos << beta(3) * (T.maxCoeff() - T.minCoeff())
              << std::noshowpos << " W\n";

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nCorrelaciones con P:  V=" << correl(V, y)
              << "  I=" << correl(Ic, y) << "  T=" << correl(T, y) << "\n";
    std::cout << "Correlaciones entre regresores: V-I=" << correl(V, Ic)
              << "  V-T=" << correl(V, T) << "  I-T=" << correl(Ic, T)
              << "  (casi nulas: sin colinealidad)\n";

    //  ENTREGABLE 5: numero de condicion de X
    std::cout << "\n===== 5. CONDICIONAMIENTO DE X =====\n";
    const VectorXd sX = svd.singularValues();
    std::cout << "Valores singulares de X: " << std::fixed << std::setprecision(6)
              << sX.transpose() << "\n";
    std::cout << "cond_2(X)   = " << cond2(X)   << "\n";
    std::cout << "cond_2(X'X) = " << cond2(XtX) << "   (= cond(X)^2 = "
              << cond2(X) * cond2(X) << ")\n";
    std::cout << "Digitos perdidos con QR           : " << std::log10(cond2(X))   << "\n";
    std::cout << "Digitos perdidos con ec. normales : " << std::log10(cond2(XtX)) << "\n";
    std::cout << "Por eso QR (o SVD) es preferible: formar X'X ELEVA AL CUADRADO\n"
                 "el numero de condicion del problema.\n";

    MatrixXd Xc(n, p);
    Xc.col(0) = VectorXd::Ones(n);
    Xc.col(1) = V.array()  - V.mean();
    Xc.col(2) = Ic.array() - Ic.mean();
    Xc.col(3) = T.array()  - T.mean();
    std::cout << "cond_2(X) centrando columnas = " << cond2(Xc)
              << "  (era " << cond2(X) << ")\n";
    std::cout << "El mal condicionamiento viene de la ESCALA dispar de las\n"
                 "columnas, no de colinealidad real entre V, I y T.\n";

    //  diagnostico fisico adicional
    VectorXd Pvi = V.array() * Ic.array();
    MatrixXd Xvi(n, 2);
    Xvi.col(0) = VectorXd::Ones(n);
    Xvi.col(1) = Pvi;
    VectorXd bvi  = Xvi.householderQr().solve(y);
    VectorXd rvi  = y - Xvi * bvi;
    const double R2vi = 1.0 - rvi.squaredNorm() / SST;

    std::cout << "\n===== DIAGNOSTICO ADICIONAL: EL MODELO FISICO =====\n";
    std::cout << "corr(P, V*I) = " << correl(y, Pvi) << "\n";
    std::cout << "Ajuste P = a + b*(V*I):  a = " << bvi(0) << ", b = " << bvi(1) << "\n";
    std::cout << "   RMSE = " << std::sqrt(rvi.squaredNorm() / n)
              << " W    R^2 = " << R2vi << "   <-- mejor que el modelo lineal\n";
    std::cout << "La potencia real es P = eta*V*I (producto, no suma). El modelo\n"
                 "lineal es una LINEALIZACION alrededor del punto de operacion medio;\n"
                 "b ~ " << bvi(1) << " es la eficiencia media del convertidor.\n";

    //  efecto del numero de muestras
    std::cout << "\n===== EFECTO DEL NUMERO DE MUESTRAS =====\n";
    std::cout << std::right << std::setw(6) << "n" << std::setw(13) << "b0"
              << std::setw(11) << "b1(V)" << std::setw(11) << "b2(I)"
              << std::setw(11) << "b3(T)" << std::setw(11) << "RMSE"
              << std::setw(11) << "cond(X)" << "\n";
    for (int m : {10, 20, 50, 100}) {
        if (m > n) break;
        MatrixXd Xm = X.topRows(m);
        VectorXd ym = y.head(m);
        VectorXd bm = Xm.householderQr().solve(ym);
        VectorXd rm = ym - Xm * bm;
        std::cout << std::setw(6) << m << std::setw(13) << bm(0) << std::setw(11) << bm(1)
                  << std::setw(11) << bm(2) << std::setw(11) << bm(3)
                  << std::setw(11) << std::sqrt(rm.squaredNorm() / m)
                  << std::setw(11) << cond2(Xm) << "\n";
    }
    std::cout << "Los coeficientes se estabilizan, los errores estandar caen como\n"
                 "1/sqrt(n) y el RMSE converge al ruido real del experimento.\n";

    return 0;
}
