#from __future__ import annotations
from pydantic import BaseModel, field_validator
from utils import unit_symbol


class ComponentBase(BaseModel):
    name: str
    value: float = 0.0
    nodes: list = []
    args: dict = {}

    @field_validator("value", mode="before")
    @classmethod
    def validate_value(cls, v: float) -> float:
        return unit_symbol(v)


class Inductor(ComponentBase):
    type_: str = "L"


class Resistor(ComponentBase):
    type_: str = "R"


class Capacitor(ComponentBase):
    type_: str = "C"


class Mos(ComponentBase):
    type_: str = "M"
    mname: str


class Current(ComponentBase):
    type_: str = "I"
    ac: float = 0.0
    dc: float = 0.0
    pulse: list[float] = [0.0] * 8
    sin: list[float] = [0.0] * 6

    @field_validator("ac", "dc", mode="before")
    @classmethod
    def unit_conversion(cls, v: float) -> float:
        return unit_symbol(v)

    @field_validator("sin", mode="before")
    @classmethod
    def sin_unit_conversion(cls, sin: list[float]) -> list[float]:
        return [unit_symbol(v) for v in sin] + [0.0] * (6 - len(sin))
    
    @field_validator("pulse", mode="before")
    @classmethod
    def pulse_unit_conversion(cls, pulse: list[float]) -> list[float]:
        return [unit_symbol(v) for v in pulse] + [0.0] * (8 - len(pulse))


class Voltage(ComponentBase):
    type_: str = "V"
    ac: float = 0.0
    dc: float = 0.0
    pulse: list[float] = [0.0] * 8
    sin: list[float] = [0.0] * 6
    pwl: list[float]=[]
    @field_validator("ac", "dc", mode="before")
    @classmethod
    def unit_conversion(cls, v: float) -> float:
        return unit_symbol(v)

    @field_validator("sin", mode="before")
    @classmethod
    def sin_unit_conversion(cls, sin: list[float]) -> list[float]:
        return [unit_symbol(v) for v in sin] + [0.0] * (6 - len(sin))
    
    @field_validator("pulse", mode="before")
    @classmethod
    def pulse_unit_conversion(cls, pulse: list[float]) -> list[float]:
        return [unit_symbol(v) for v in pulse] + [0.0] * (8 - len(pulse))
        #V1 V2 TD TR TF PW PER NP
    def pwl_unit_conversion(cls, pulse: list[float]) -> list[float]:
        
        return [unit_symbol(v) for v in pulse] + [0.0] * (8 - len(pulse))

class Diode(ComponentBase):
    type_: str = "D"
    mname: str
