// Унифицированный интерфейс дифференциального уравнения
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    /// <summary>
    /// Унифицированный интерфейс дифференциального уравнения
    /// </summary>
    public interface IDifferentialEquation
    {
        /// <summary>
        /// Крайняя левая точка диапазона "х" 
        /// </summary>
        double Min_X { get; set; }
        /// <summary>
        /// Значение "у" в крайней левой точке 
        /// </summary>
        double  Min_Y { get; set; }
        /// <summary>
        /// Значение "z" в крайней левой точке 
        /// </summary>
        double Min_Z { get; set; }
        /// <summary>
        /// Шаг сетки 
        /// </summary>
        double Step { get; set; }
        /// <summary>
        /// Количество итераций c шагом Step  
        /// </summary>
        int CountIterations { get; set; }
        /// <summary>
        /// Округление до нужного знака, после запятой
        /// </summary>
        int Rounding { get; set; }
        /// <summary>
        /// Рассчитать хранимое уравнение
        /// </summary>
        /// <param name="Values">Значения</param>
        /// <returns>Результат</returns>
        double ComputeEquation(params double[] Values);
    }
}


