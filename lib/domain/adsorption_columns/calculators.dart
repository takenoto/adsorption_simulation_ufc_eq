import 'dart:math' as math;

///Representa a variação da concentração com o tempo
class UnidimensionalDcDt {
  ///U é a velocidade linear intersticial
  double dc_dt(
      {required double dax,
      required double d2c_dx2,
      required double nu,
      required double u,
      required double dc_dx,
      required double porosity,
      required double dq_dt}) {
    return dax * d2c_dx2 - u * dc_dx + ((1 - porosity) / porosity) * dq_dt;
  }
}

typedef dq_dt_FromParamsBase = double Function({required double C, required double q});
class LumpedDqDtCalculator<T>{
  //[T] é a key, que representa a chave para a substância desse adsorvente que tem esse comportamento.
  Map<T, dq_dt_FromParamsBase> _map;
  
  dq_dt_FromParamsBase forSubstance(T key) => _map[key]!;

  LumpedDqDtCalculator(this._map);

  
}

class ValuePair {
  double v1;
  double v2;
  ValuePair(this.v1, this.v2);
}

///Derivada do tipo dy/dx
class SingleDerivativeCalculatorFiniteDifference {
  //Calcula derivada simples com os resultados de y e x passados
  //Pode ser usado tanto para diferença forward qunanto backward. Até pra central também!
  double difference(ValuePair y, ValuePair x) {
    return (y.v2 - y.v1) / (x.v2 - x.v1);
  }
}

class SecondDerivativeInput {
  ///O valor da função após o ponto que está sendo avaliado
  double f_plus_1;

  ///O valor da função no ponto em que a derivada é avaliada
  double f_0;

  ///O valor da função no ponto anterior ao em que a derivada é avaliada
  double f_minus_1;

  ///A diferença entre a variável que deriva e o próximo valor e o anterior.
  ///
  ///Ex: (f(x + h) - 2*f(x) + f(x - h))/h²
  double h;
  SecondDerivativeInput(
      {required this.f_minus_1,
      required this.f_0,
      required this.f_plus_1,
      required this.h});
}

class SecondDerivativeCalculatorFiniteDifference {
  ///Cálculo para derivada do tipo d²y/dx²
  double d2ydx2(SecondDerivativeInput i) {
    return (i.f_plus_1 - 2 * i.f_0 + i.f_minus_1) / (math.pow(i.h, 2));
  }
}
