// Делегат дифференциальных уравнений вида dy/dx = f(x, y)
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Делегат дифференциального уравнения вида dy/dx = f(x, y)
    /// </summary>
    /// <param name="X">Значение X</param>
    /// <param name="Y">Значение Y</param>
    /// <returns>Значение производной</returns>
    public delegate double AEquation_dYdX(double X, double Y);
 
}
