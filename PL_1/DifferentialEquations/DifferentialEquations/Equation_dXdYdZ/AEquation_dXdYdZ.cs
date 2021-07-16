//Делегат дифференциального уравнения вида d^2y/dx^2+p*dy/dx + q = f(x, y) 
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Делегат дифференциального уравнения вида d^2y/dx^2+p*dy/dx + q = f(x,y) 
    /// </summary>
    /// <param name="X">Значение X</param>
    /// <param name="Y">Значение Y</param>
    /// <param name="Z">Значение Z</param>
    /// <returns>Значение производной</returns>
    public delegate double AEquation_dXdYdZ(double X, double Y, double Z);
   
}