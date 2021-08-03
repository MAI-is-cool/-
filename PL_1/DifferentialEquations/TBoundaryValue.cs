//Класс описывающий свойства граничного значения, по которому будет прекращаться интегрирование.
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
//*********************************************************
namespace StandartHelperLibrary.MathHelper
{
    public class TBoundaryValue
    {
//--------------------------------------------------------------------------------------------
        /// <summary>
        /// Название граничной величины
        /// </summary>
        public string Name { get; set; }
        /// <summary>
        /// Значение граничной величины
        /// </summary>
        public double BoundaryValue { get; set; }
        /// <summary>
        /// Потребная точность невязки на выходе
        /// </summary>
        public double Accuracy { get; set; }
        /// <summary>
        /// Текущее значение граничной величены
        /// </summary>
        public double CurrentValue { get; set; }
        /// <summary>
        /// Возвращате значение невязки на текущей итерации
        /// </summary>
        public double Residual
        {
            get
            {
                return BoundaryValue - CurrentValue;
            }
                
        }
//--------------------------------------------------------------------------------------------
    }
}
