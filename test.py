from flask import Flask, render_template, request
from scipy.optimize import linprog
import numpy as np

app = Flask(__name__)


def data_const():
    data = {}
    #  data[''] = 0
    #  data['A1'] = 0
    data['sconst_f'] = 0.847

    # �������� �� ���������� ��������� (� ����)
    data['mn1'] = 19
    data['mn6'] = 36.6
    data['mn7'] = 60
    data['mn8'] = 9
    data['mn9'] = 1.88
    data['mn10'] = 10
    data['mn11'] = 10

    # ��������� �� ���������� ��������� (� ����)
    data['mx1'] = 160
    data['mx6'] = 36.9
    data['mx7'] = 80
    data['mx8'] = 19.5
    data['mx9'] = 5.7
    data['mx10'] = 60
    data['mx11'] = 140

    data['sx1'] = 0
    data['sx2'] = 0
    data['sx3'] = 0
    data['sx4'] = 0
    data['sx5'] = 0
    data['sx6'] = 0
    data['sx7'] = 0
    data['sx8'] = 0
    data['sx9'] = 0
    data['sx10'] = -0.0310  # ????
    data['sx11'] = 0.033  # ?????

    #######################################
    #    data['   ']   =

    data['A1'] = -0.002
    data['A2'] = -0.210
    data['A3'] = 0.001
    data['A4'] = 0.003
    data['A5'] = -0.057
    data['A6'] = 0
    data['A7'] = 0
    data['A8'] = 0
    data['A9'] = 0
    data['A10'] = 0.002
    data['A11'] = 0.000

    data['A21'] = -0.029
    data['A22'] = 9.348
    data['A23'] = 0.271
    data['A24'] = -0.498
    data['A25'] = 2.178
    data['A26'] = 0
    data['A27'] = 0
    data['A28'] = 0
    data['A29'] = 0
    data['A210'] = 0.109
    data['A211'] = -0.118

    data['A31'] = -0.022
    data['A32'] = 1.569
    data['A33'] = -0.025
    data['A34'] = -0.013
    data['A35'] = -0.741
    data['A36'] = 0
    data['A37'] = 0
    data['A38'] = 0
    data['A39'] = 0
    data['A310'] = 0.025
    data['A311'] = -0.009

    data['A1const'] = 44.885
    data['A2const'] = -292.464
    data['A3const'] = -31.555
    return data


def linpr(pacient, data=data_const()):
    # Standart values button

    const_f = data['sconst_f']

    edit_age = pacient["age"]
    edit22 = pacient["temperature"]  # Temperatura
    edit23 = pacient["edit23"]  # CHSS
    edit24 = pacient["edit24"]  # Tiroksin
    edit25 = pacient["edit25"]  # Triyuoditron

    data['f'] = np.asarray(
        [data['sx1'], data['sx2'], data['sx3'], data['sx4'], data['sx5'], data['sx6'], data['sx7'], data['sx8'],
         data['sx9'], data['sx10'], data['sx11']])
    data['A'] = [[data['A1'], data['A2'], data['A3'], data['A4'], data['A5'], 0, 0, 0, 0, data['A10'], data['A11']],
                 [data['A21'], data['A22'], data['A23'], data['A24'], data['A25'], 0, 0, 0, 0, data['A210'],
                  data['A211']],
                 [data['A31'], data['A32'], data['A33'], data['A34'], data['A35'], 0, 0, 0, 0, data['A310'],
                  data['A311']],
                 [- data['A1'], - data['A2'], - data['A3'], - data['A4'], - data['A5'], 0, 0, 0, 0, - data['A10'],
                  - data['A11']],
                 [- data['A21'], - data['A22'], - data['A23'], - data['A24'], - data['A25'], 0, 0, 0, 0, - data['A210'],
                  - data['A211']],
                 [- data['A31'], - data['A32'], - data['A33'], - data['A34'], - data['A35'], 0, 0, 0, 0, - data['A310'],
                  - data['A311']]]

    data['intcon'] = 1
    data['b'] = [37 - 44.885, 80 + 292.464, 19.5 + 31.555, -36 + 44.885, -50 - 292.464, -9 - 31.555]
    data['ub'] = np.asarray(
        [data['mx1'], 43, 220, 40, 10, data['mx6'], data['mx7'], data['mx8'], data['mx9'], data['mx10'], data['mx11']])
    data['lb'] = np.asarray(
        [data['mn1'], 0, 0, 0, 0, data['mn6'], data['mn7'], data['mn8'], data['mn9'], data['mn10'], data['mn11']])

    # data['ub'] = [[data['mx1']], [43], [220], [40], [10], [data['mx6']] ,[data['mx7']], [data['mx8']], [data['mx9']],[data['mx10']], [data['mx11']]]
    # data['lb'] = [[data['mn1']], [0], [0], [0], [0], [data['mn6']], [data['mn7']], [data['mn8']], [data['mn9']], [data['mn10']], [data['mn11']]]

    data['Aeq'] = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                   [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]]

    A = np.array([1, 2]).reshape([1, 2])
    zerosm = np.zeros([1, 2])
    A_eq = np.hstack([A, zerosm])

    data['beq'] = [edit_age, edit22, edit23, edit24,
                   edit25]  # data['edit_age data['edit22 data['edit23 data['edit24 data['edit25

    res = linprog(c=data['f'], A_ub=data['A'], b_ub=data['b'], bounds=list(zip(data['lb'], data['ub'])),
                  A_eq=data['Aeq'], b_eq=data['beq'], options={"disp": True})
    # , data['beq']), data['lb'],data['ub'])   bounds= (data['lb'],data['ub'])
    # res = linprog(data['f'],data['intcon'],data['A'], data['b'], data['Aeq'], data['beq'], data['lb'],data['ub'])

    answer = res['fun']
    dr = res['x']
    dr1 = dr[-2]
    dr2 = dr[-1]
    print(dr1, dr2, answer)
    return {"dr1": dr1, "dr2": dr2, "answer": answer, "data": data}


@app.route("/", methods=["POST", "GET"])
def index():
    return render_template("page2.html")


@app.route("/main", methods=["POST", "GET"])
def Main():
    pacient = {"age": 0, "temperature": 0, "edit23": 0, "edit24": 0, "edit25": 0}
    count = 0
    for k in pacient.keys():
        if request.args.get(k):
            pacient[k] = float(request.args.get(k))
            count += 1

    if count < 3:
        return render_template("index.html", data=data_const(), pacient=pacient)

    if request.args.get("ok2"):
        lp = linpr(pacient)
        data = lp['data']
        return render_template("index.html", data=data, lp=lp, pacient=pacient)

        # 1
    dima = data_const()
    for s in dima.keys():
        if request.args.get(s):
            dima[s] = float(request.args.get(s))
    lp = linpr(pacient, dima)
    return render_template("index.html", lp=lp, data=dima, pacient=pacient)


if __name__ == "__main__":
    app.run(debug=True)

