#pragma once

#include <TMatrixD.h>
#include <TLorentzVector.h>

#include <utility>


/**
 * \class NuRecoRochester
 * \brief This class reconstructs neutrino from a t -> blv decay using Rochester algorithm
 * 
 * The algorithm is described in [1]. It applies constraints from masses of the top quark and the
 * W boson, which lead to an ellipsis in the space of neutrino three-momentum. A unique solution on
 * this ellipsis is chosen based on the compatibility with the measured MET.
 * [1] B.A. Betchart, R. Demina, A. Harel, Nucl.Instrum.Meth. A736 (2014) 169 [arXiv:1305.1878]
 * 
 * This method was used in TOP-16-008 (AN-16-020).
 * 
 * Implementation is copied from [2-3].
 * [2] https://gitlab.cern.ch/mverzett/URTTbar/blob/fd3362a007bdc0bea3f9136dff6de1a700645488/interface/NeutrinoSolver.h
 * [3] https://gitlab.cern.ch/mverzett/URTTbar/blob/fd3362a007bdc0bea3f9136dff6de1a700645488/src/NeutrinoSolver.cc
 */
class NuRecoRochester
{
public:
    /**
     * \brief Constructor from a lepton and a b-quark jet
     * 
     * If netrino cannot be reconstructed for this configuration (which, presumably, happens when
     * m(l, b) > mt), an error flag is set.
     */
    NuRecoRochester(TLorentzVector const *lep, TLorentzVector const *bjet,
      double MW = 80, double MT = 173);
    
public:
    /**
     * \brief Finds neutrino solution that minimizes figure of merit computed by method Chi2
     * 
     * Parameter test will contain the computed figure of merit. Version of the algorithm described
     * in the paper is reproduced with an identity MET error matrix (metxerr = metyerr = 1,
     * metxyrho = 0.).
     * 
     * If solution cannot be found, returns a zero four-vector and test is set to -1.
     */
    TLorentzVector GetBest(double metx, double mety, double metxerr, double metyerr,
      double metxyrho, double &test, bool INFO = false);
    
private:
    /// Constructs a matrix for rotation about x axis through an angle a
    TMatrixD RotationX(double a);
    
    /// Constructs a matrix for rotation about y axis through an angle a
    TMatrixD RotationY(double a);
    
    /// Constructs a matrix for rotation about z axis through an angle a
    TMatrixD RotationZ(double a);
    
    /**
     * \brief Constructs neutrino solution for the given parameter t
     * 
     * The argument seems to be a parameter that identifies a point on the solution ellipsis.
     */
    void Solve(double t);
    
    /**
     * \brief Constructs neutrino solution and reports it as a vector in the transverse plane
     * 
     * The solution is returned in the form of a 2x1 matrix. Consult documentation for method
     * Solve for the meaning of the parameter t.
     */
    TMatrixD GetPtSolution(double t);
    
    /**
     * \brief Constructs neutrino solution and reports it as a four-vector
     * 
     * Consult documentation for method Solve for the meaning of the parameter t.
     */
    TLorentzVector GetSolution(double t);
    
    /**
     * \Computes figure of merit for the neutrino solution
     * 
     * Consult documentation for method Solve for the meaning of parameter t. When the MET error
     * matrix is identity, the returned value is the Euclidian distance in the transverse plane
     * between the neutrino solution given by the parameter t and experimental MET.
     */
    double Chi2(double t);
    
    /**
     * \brief Finds extremum of function Chi2
     * 
     * Finds minimum if MIN is true and maximum otherwise. Returns a pair consisting of a point of
     * the extremum and value of Chi2 evaluated there.
     */
    std::pair<double, double> Extrem(double t, bool MIN = true);
    
private:
    /// Masses of top quark, W boson, lepton, b-quark jet, and neutrino (Mn = 0.)
    double Mt, Mw, Ml, Mb, Mn;
    
    /// Error flag set when no solution can be found for given b-quark jet and lepton
    bool ERROR;
    
    TMatrixD H;
    
    /// Current solution for neutrino three-momentum (3x1 matrix)
    TMatrixD T;
    
    /// Experimental MET stored as a 2x1 matrix
    TMatrixD MET;
    
    /// Iverted MET error matrix (2x2)
    TMatrixD VM;
};
