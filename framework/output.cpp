#include <cmath>
#include <Eigen/Dense>

double distanceNorm(Eigen::Vector3d x1, Eigen::Vector3d x2) {
    return (x1-x2).norm();
}

double s(double q) {
    return ((3)/(2*(M_PI)))*
    ((0 <= q && q < 1) ? ((2)/(3))-std::pow(q, 2)+((1)/(2))*std::pow(q, 3) : ((1 <= q && q < 2) ? ((1)/(6))*(2-q)^3 : 0));
}
double W(Eigen::Vector3d x1, Eigen::Vector3d x2, double h) {
    return ((1)/(std::pow(h, 3)))*s((((distanceNorm(x1,x2)))/(h)));
}
double sgrad(double q) {
    return ((0 <= q && q < 1) ? (((3)/(2*(M_PI)))) * (-2*q+((3)/(2))*std::pow(q, 2)) : ((1 <= q && q < 2) ? (((3)/(2*(M_PI)))) * (((-1)/(2))*(2-q)^2) : 0));
}
double Wgrad(Eigen::Vector3d x1, Eigen::Vector3d x2, double h) {
    return ((x1 == x2) ? { 0 ,  0 ,  0 } : {
            ((((1)/(std::pow(h, 4)))*sgrad((((distanceNorm(x1,x2)))/(h))))/((distanceNorm(x1,x2)))) * (x1.x - x2.x) , 
            ((((1)/(std::pow(h, 4)))*sgrad((((distanceNorm(x1,x2)))/(h))))/((distanceNorm(x1,x2)))) * (x1.y - x2.y) , 
            ((((1)/(std::pow(h, 4)))*sgrad((((distanceNorm(x1,x2)))/(h))))/((distanceNorm(x1,x2)))) * (x1.z - x2.z)
        });
}
double xsph(double m, Eigen::Vector3d x1, Eigen::Vector3d x2, Eigen::Vector3d v1, Eigen::Vector3d v2, double p1, double p2, double h) {
    return 2*m*((v2-v1)/(p1+p2))*W(x1,x2,h);
}
double i(double B, double pi, double p0) {
    return max(0, B*(pi-p0));
}
double apfluidtofluid(double m, Eigen::Vector3d x1, Eigen::Vector3d x2, double pressure1, double pressure2, double density1, double density2, double h) {
    return m*(((pressure1)/(std::pow(density1, 2)))+((pressure2)/(std::pow(density2, 2))))*Wgrad(x1,x2,h);
}
double apfluidtoboundary(double restDens, Eigen::Vector3d x1, Eigen::Vector3d x2, double volume, double pressure1, double density1, double h) {
    return restDens*volume*(((pressure1)/(std::pow(density1, 2))))*Wgrad(x1,x2,h);
}
double aviscosfluidtofluid(Eigen::Vector3d x1, Eigen::Vector3d x2, Eigen::Vector3d v1, Eigen::Vector3d v2, double m, double density2, double h) {
    return ((m)/(density2))*(v1-v2)*((((x1-x2)*Wgrad(x1,x2,h))/((distanceNorm(x1,x2))*(distanceNorm(x1,x2))+0.01*std::pow(h, 2))));
}
double aviscosfluidtoboundary(Eigen::Vector3d x1, Eigen::Vector3d x2, Eigen::Vector3d v1, double volume, double h) {
    return volume*v1*((((x1-x2)*Wgrad(x1,x2,h))/((distanceNorm(x1,x2))*(distanceNorm(x1,x2))+0.01*std::pow(h, 2))));
}
double Wcoh(double r, double c) {
    return ((32)/((M_PI)*std::pow(c, 9)))*
    ((0 <= r && r <= ((c)/(2))) ? 2*(c-r)^3*std::pow(r, 3)-((std::pow(c, 6))/(64)) : ((((c)/(2)) < r && r <= c) ? (c-r)^3*std::pow(r, 3) : 0));
}
double Wadh(double r, double c) {
    return ((0.007)/(std::pow(c, 3.25))*(((c)/(2) <= r && r <= c) ? (((-4*std::pow(r, 2))/(c))+6*r-2*c)^{1/4} : 0);
}
double fcohesion(double cohFac, double m, Eigen::Vector3d x1, Eigen::Vector3d x2, double c) {
    return -cohFac*m*m*Wcoh((distanceNorm(x1,x2)),c)*((x1-x2)/((distanceNorm(x1,x2))));
}
double ni(double c, double m, double density, Eigen::Vector3d x1, Eigen::Vector3d x2, double h) {
    return c*((m)/(density))*Wgrad(x1,x2,h);
}
double fcurvature(double cohFac, double m, Eigen::Vector3d n1, Eigen::Vector3d n2) {
    return -cohFac*m*(n1-n2);
}
double fadhesion(double adhFac, double mFluid, double mBound, Eigen::Vector3d x1, Eigen::Vector3d x2, double c) {
    return -adhFac*mFluid*mBound*Wadh((distanceNorm(x1,x2)),c)*((x1-x2)/((distanceNorm(x1,x2))));
}
