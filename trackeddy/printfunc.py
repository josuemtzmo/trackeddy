import sys
import time


class Printer:
    """
    Print things to stdout on one line dynamically
    """

    def __init__(self):
        import platform

        if platform.python_version_tuple() < ("3", "10", "0"):
            raise ImportError("To use this option update to python > 3.10")
        self.tic = time.time()
        sys.stdout.flush()
        self.data = []

    def printtextoneline(self, string):
        sys.stdout.write("\r\x1b[K" + string.__str__())
        sys.stdout.flush()

    def timepercentprint(
        self,
        minv,
        maxv,
        step,
        i,
        neddies,
    ):
        percent = (float(i + 1) / maxv) * 100.0
        etime = round(time.time() - self.tic)
        stmtime = round((etime / percent) * 100)

        progress = int(10 / (maxv / (step * (i + 1))))
        emptyprog = 10 - progress

        sys.stdout.write(
            """\r 0% [{0}>{1}]{2}% | Elapsed Time: {3} s |
            Estimated Time: {4} s | Info: {5} |""".format(
                "=" * progress,
                " " * emptyprog,
                round(percent),
                etime,
                stmtime,
                neddies,
            )
        )
        if percent != 100 and i != maxv:
            sys.stdout.flush()
        else:
            # print(self.data)
            # plt.plot(self.data)
            # plt.show()
            print("")
