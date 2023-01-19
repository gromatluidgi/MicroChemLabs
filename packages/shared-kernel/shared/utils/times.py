from datetime import datetime, timedelta


class DatetimeUtils:
    @staticmethod
    def yesterday():
        return datetime.now() - timedelta(days=1)
