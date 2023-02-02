function Euler2Quaternion = Euler2Quaternion(phi, theta, psi)

Euler2Quaternion(1) =  cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2);

Euler2Quaternion(2) =  cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2);

Euler2Quaternion(3) =  cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2);

Euler2Quaternion(4) =  sin(psi/2)*cos(theta/2)*cos(phi/2) - cos(psi/2)*sin(theta/2)*sin(phi/2);